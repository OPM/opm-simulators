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
 * \brief Routines that actually solve the ODEs that emerge from the hydrostatic
 *        equilibrium problem
 */
#ifndef EWOMS_INITSTATEEQUIL_HH
#define EWOMS_INITSTATEEQUIL_HH

#include <opm/models/utils/propertysystem.hh>

#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/fluidstates/SimpleModularFluidState.hpp>

#include <array>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>
#include <string>

namespace Opm {

class EclipseState;
class EquilRecord;
class NumericalAquifers;

/**
 * Types and routines that collectively implement a basic
 * ECLIPSE-style equilibration-based initialisation scheme.
 *
 * This namespace is intentionally nested to avoid name clashes
 * with other parts of OPM.
 */
namespace EQUIL {

class EquilReg;
namespace Miscibility { class RsFunction; }

namespace Details {
template <class RHS>
class RK4IVP {
public:
    RK4IVP(const RHS& f,
           const std::array<double,2>& span,
           const double y0,
           const int N);

    double
    operator()(const double x) const;

private:
    int N_;
    std::array<double,2> span_;
    std::vector<double>  y_;
    std::vector<double>  f_;

    double stepsize() const;
};

namespace PhasePressODE {
template <class FluidSystem>
class Water
{
using TabulatedFunction = Tabulated1DFunction<double>;
public:
    Water(const TabulatedFunction& tempVdTable,
          const TabulatedFunction& saltVdTable,
          const int pvtRegionIdx,
          const double normGrav);

    double operator()(const double depth,
                      const double press) const;

private:
    const TabulatedFunction& tempVdTable_;
    const TabulatedFunction& saltVdTable_;
    const int pvtRegionIdx_;
    const double g_;

    double density(const double depth,
                   const double press) const;
};

template <class FluidSystem, class RS>
class Oil
{
using TabulatedFunction = Tabulated1DFunction<double>;
public:
    Oil(const TabulatedFunction& tempVdTable,
        const RS& rs,
        const int pvtRegionIdx,
        const double normGrav);

    double operator()(const double depth,
                      const double press) const;

private:
    const TabulatedFunction& tempVdTable_;
    const RS& rs_;
    const int pvtRegionIdx_;
    const double g_;

    double density(const double depth,
                   const double press) const;
};

template <class FluidSystem, class RV, class RVW>
class Gas
{
using TabulatedFunction = Tabulated1DFunction<double>;
public:
    Gas(const TabulatedFunction& tempVdTable,
        const RV& rv,
        const RVW& rvw,
        const int pvtRegionIdx,
        const double normGrav);

    double operator()(const double depth,
                      const double press) const;

private:
    const TabulatedFunction& tempVdTable_;
    const RV& rv_;
    const RVW& rvw_;
    const int pvtRegionIdx_;
    const double g_;

    double density(const double depth,
                   const double press) const;
};

} // namespace PhasePressODE

template <class FluidSystem, class Region>
class PressureTable
{
public:
    using VSpan = std::array<double, 2>;

    /// Constructor
    ///
    /// \param[in] gravity Norm of gravity vector (acceleration strength due
    ///    to gravity).  Normally the standardised value at Tellus equator
    ///    (9.80665 m/s^2).
    ///
    /// \param[in] samplePoints Number of equally spaced depth sample points
    ///    in each internal phase pressure table.
    explicit PressureTable(const double gravity,
                           const int    samplePoints = 2000);

    /// Copy constructor
    ///
    /// \param[in] rhs Source object for copy initialization.
    PressureTable(const PressureTable& rhs);

    /// Move constructor
    ///
    /// \param[in,out] rhs Source object for move initialization.  On output,
    ///    left in a moved-from ("valid but unspecified") state.  Internal
    ///    pointers in \p rhs are null (\c unique_ptr guarantee).
    PressureTable(PressureTable&& rhs);

    /// Assignment operator
    ///
    /// \param[in] rhs Source object.
    ///
    /// \return \code *this \endcode.
    PressureTable& operator=(const PressureTable& rhs);

    /// Move-assignment operator
    ///
    /// \param[in] rhs Source object.  On output, left in a moved-from ("valid
    ///    but unspecified") state.  Internal pointers in \p rhs are null (\c
    ///    unique_ptr guarantee).
    ///
    /// \return \code *this \endcode.
    PressureTable& operator=(PressureTable&& rhs);

    void equilibrate(const Region& reg,
                     const VSpan&  span);

    /// Predicate for whether or not oil is an active phase
    bool oilActive() const;

    /// Predicate for whether or not gas is an active phase
    bool gasActive() const;

    /// Predicate for whether or not water is an active phase
    bool waterActive() const;

    /// Evaluate oil phase pressure at specified depth.
    ///
    /// \param[in] depth Depth of evaluation point.  Should generally be
    ///    within the \c span from the previous call to \code equilibrate()
    ///    \endcode.
    ///
    /// \return Oil phase pressure at specified depth.
    double oil(const double depth) const;

    /// Evaluate gas phase pressure at specified depth.
    ///
    /// \param[in] depth Depth of evaluation point.  Should generally be
    ///    within the \c span from the previous call to \code equilibrate()
    ///    \endcode.
    ///
    /// \return Gas phase pressure at specified depth.
    double gas(const double depth) const;

    /// Evaluate water phase pressure at specified depth.
    ///
    /// \param[in] depth Depth of evaluation point.  Should generally be
    ///    within the \c span from the previous call to \code equilibrate()
    ///    \endcode.
    ///
    /// \return Water phase pressure at specified depth.
    double water(const double depth) const;

private:
    template <class ODE>
    class PressureFunction
    {
    public:
        struct InitCond {
            double depth;
            double pressure;
        };

        explicit PressureFunction(const ODE&      ode,
                                  const InitCond& ic,
                                  const int       nsample,
                                  const VSpan&    span);

        PressureFunction(const PressureFunction& rhs);

        PressureFunction(PressureFunction&& rhs) = default;

        PressureFunction& operator=(const PressureFunction& rhs);

        PressureFunction& operator=(PressureFunction&& rhs);

        double value(const double depth) const;

    private:
        enum Direction : std::size_t { Up, Down, NumDir };

        using Distribution = Details::RK4IVP<ODE>;
        using DistrPtr = std::unique_ptr<Distribution>;

        InitCond initial_;
        std::array<DistrPtr, Direction::NumDir> value_;
    };

    using OilPressODE = PhasePressODE::Oil<
        FluidSystem, typename Region::CalcDissolution
    >;

    using GasPressODE = PhasePressODE::Gas<
        FluidSystem, typename Region::CalcEvaporation, typename Region::CalcWaterEvaporation
    >;

    using WatPressODE = PhasePressODE::Water<FluidSystem>;

    using OPress = PressureFunction<OilPressODE>;
    using GPress = PressureFunction<GasPressODE>;
    using WPress = PressureFunction<WatPressODE>;

    using Strategy = void (PressureTable::*)
        (const Region&, const VSpan&);

    double gravity_;
    int    nsample_;

    std::unique_ptr<OPress> oil_{};
    std::unique_ptr<GPress> gas_{};
    std::unique_ptr<WPress> wat_{};

    template <typename PressFunc>
    void checkPtr(const PressFunc*   phasePress,
                  const std::string& phaseName) const;

    Strategy selectEquilibrationStrategy(const Region& reg) const;

    void copyInPointers(const PressureTable& rhs);

    void equil_WOG(const Region& reg, const VSpan& span);
    void equil_GOW(const Region& reg, const VSpan& span);
    void equil_OWG(const Region& reg, const VSpan& span);

    void makeOilPressure(const typename OPress::InitCond& ic,
                         const Region&                    reg,
                         const VSpan&                     span);

    void makeGasPressure(const typename GPress::InitCond& ic,
                         const Region&                    reg,
                         const VSpan&                     span);

    void makeWatPressure(const typename WPress::InitCond& ic,
                         const Region&                    reg,
                         const VSpan&                     span);
};

// ===========================================================================

/// Simple set of per-phase (named by primary component) quantities.
struct PhaseQuantityValue {
    double oil{0.0};
    double gas{0.0};
    double water{0.0};

    PhaseQuantityValue& axpy(const PhaseQuantityValue& rhs, const double a)
    {
        this->oil   += a * rhs.oil;
        this->gas   += a * rhs.gas;
        this->water += a * rhs.water;

        return *this;
    }

    PhaseQuantityValue& operator/=(const double x)
    {
        this->oil   /= x;
        this->gas   /= x;
        this->water /= x;

        return *this;
    }

    void reset()
    {
        this->oil = this->gas = this->water = 0.0;
    }
};

/// Calculator for phase saturations
///
/// Computes saturation values at arbitrary depths.
///
/// \tparam MaterialLawManager Container for material laws.  Typically a
///    specialization of the \code Opm::EclMaterialLawManager<> \endcode
///    template.
///
/// \tparam FluidSystem An OPM fluid system type.  Typically a
///    specialization of the \code Opm::BlackOilFluidSystem<> \endcode
///    template.
///
/// \tparam Region Representation of an equilibration region.  Typically
///    \code Opm::EQUIL::EquilReg \endcode from the equilibrationhelpers.
///
/// \tparam CellID Representation an equilibration region's cell IDs.
///    Typically \code std::size_t \endcode.
template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
class PhaseSaturations
{
public:
    /// Evaluation point within a model geometry.
    ///
    /// Associates a particular depth to specific cell.
    struct Position {
        CellID cell;
        double depth;
    };

    /// Convenience type alias
    using PTable = PressureTable<FluidSystem, Region>;

    /// Constructor
    ///
    /// \param[in,out] matLawMgr Read/write reference to a material law
    ///    container.  Mutated by member functions.
    ///
    /// \param[in] swatInit Initial water saturation array (from SWATINIT
    ///    data).  Empty if SWATINIT is not used in this simulation model.
    explicit PhaseSaturations(MaterialLawManager&        matLawMgr,
                              const std::vector<double>& swatInit);

    /// Copy constructor.
    ///
    /// \param[in] rhs Source object.
    PhaseSaturations(const PhaseSaturations& rhs);

    /// Disabled assignment operator.
    PhaseSaturations& operator=(const PhaseSaturations&) = delete;

    /// Disabled move-assignment operator.
    PhaseSaturations& operator=(PhaseSaturations&&) = delete;

    /// Calculate phase saturations at particular point of the simulation
    /// model geometry.
    ///
    /// \param[in] x Specific geometric point (depth within a specific cell).
    ///
    /// \param[in] reg Equilibration information for a single equilibration
    ///    region; notably contact depths.
    ///
    /// \param[in] ptable Previously equilibrated phase pressure table
    ///    pertaining to the equilibration region \p reg.
    ///
    /// \return Set of phase saturation values defined at particular point.
    const PhaseQuantityValue&
    deriveSaturations(const Position& x,
                      const Region&   reg,
                      const PTable&   ptable);

    /// Retrieve saturation-corrected phase pressures
    ///
    /// Values associated with evaluation point of previous call to \code
    /// deriveSaturations() \endcode.
    const PhaseQuantityValue& correctedPhasePressures() const
    {
        return this->press_;
    }

private:
    /// Convenience amalgamation of the deriveSaturations() input state.
    /// These values are almost always used in concert.
    struct EvaluationPoint {
        const Position* position{nullptr};
        const Region*   region  {nullptr};
        const PTable*   ptable  {nullptr};
    };

    /// Simplified fluid state object that contains only the pieces of
    /// information needed to calculate the capillary pressure values from
    /// the current set of material laws.
    using FluidState = ::Opm::
        SimpleModularFluidState<double, /*numPhases=*/3, /*numComponents=*/3,
                                FluidSystem,
                                /*storePressure=*/false,
                                /*storeTemperature=*/false,
                                /*storeComposition=*/false,
                                /*storeFugacity=*/false,
                                /*storeSaturation=*/true,
                                /*storeDensity=*/false,
                                /*storeViscosity=*/false,
                                /*storeEnthalpy=*/false>;

    /// Convenience type alias.
    using MaterialLaw = typename MaterialLawManager::MaterialLaw;

    /// Fluid system's representation of phase indices.
    using PhaseIdx = std::remove_cv_t<
        std::remove_reference_t<decltype(FluidSystem::oilPhaseIdx)>
    >;

    /// Read/write reference to client's material law container.
    MaterialLawManager& matLawMgr_;

    /// Client's SWATINIT data.
    const std::vector<double>& swatInit_;

    /// Evaluated phase saturations.
    PhaseQuantityValue sat_;

    /// Saturation-corrected phase pressure values.
    PhaseQuantityValue press_;

    /// Current evaluation point.
    EvaluationPoint evalPt_;

    /// Capillary pressure fluid state.
    FluidState fluidState_;

    /// Evaluated capillary pressures from current set of material laws.
    std::array<double, FluidSystem::numPhases> matLawCapPress_;

    /// Capture the input evaluation point information in internal state.
    ///
    /// \param[in] x Specific geometric point (depth within a specific cell).
    ///
    /// \param[in] reg Equilibration information for a single equilibration
    ///    region; notably contact depths.
    ///
    /// \param[in] ptable Previously equilibrated phase pressure table
    ///    pertaining to the equilibration region \p reg.
    void setEvaluationPoint(const Position& x,
                            const Region&   reg,
                            const PTable&   ptable);

    /// Initialize phase saturation and phase pressure values.
    ///
    /// Looks up phase pressure values from the input pressure table.
    void initializePhaseQuantities();

    /// Derive phase saturation for oil.
    ///
    /// Calculated as 1 - Sw - Sg.
    void deriveOilSat();

    /// Derive phase saturation for gas.
    ///
    /// Inverts capillary pressure curve if non-constant or uses a simple
    /// depth consideration with respect to G/O contact depth otherwise.
    void deriveGasSat();

    /// Derive phase saturation for water.
    ///
    /// Uses input data if simulation model is defined in terms of SWATINIT.
    /// Otherwise, inverts capillary pressure curve if non-constant or uses
    /// a simple depth consideration with respect to the O/W contact depth
    /// if capillary pressure curve is constant within the current cell.
    void deriveWaterSat();

    /// Correct phase saturation and pressure values to account for
    /// overlapping transition zones between G/O and O/W systems.
    void fixUnphysicalTransition();

    /// Re-adjust phase pressure values to account for phase saturations
    /// outside permissible ranges.
    void accountForScaledSaturations();

    // --------------------------------------------------------------------
    // Note: Function 'applySwatInit' is non-const because the overload set
    // needs to mutate the 'matLawMgr_'.
    // --------------------------------------------------------------------

    /// Derive water saturation from SWATINIT data.
    ///
    /// Uses SWATINIT array data from current cell directly.  Also updates
    /// the material law container's internal notion of the maximum
    /// attainable O/W capillary pressure value.
    ///
    /// \param[in] pcow O/W capillary pressure value (Po - Pw).
    ///
    /// \return Water saturation value.
    std::pair<double, bool> applySwatInit(const double pcow);

    /// Derive water saturation from SWATINIT data.
    ///
    /// Uses explicitly passed-in saturation value.  Also updates the
    /// material law container's internal notion of the maximum attainable
    /// O/W capillary pressure value.
    ///
    /// \param[in] pc x/W capillary pressure value (Px - Pw; x in {O, G}).
    ///
    /// \param[in] sw Water saturation value.
    ///
    /// \return Water saturation value.  Input value, possibly mollified by
    ///    current set of material laws.
    std::pair<double, bool> applySwatInit(const double pc, const double sw);

    /// Invoke material law container's capillary pressure calculator on
    /// current fluid state.
    void computeMaterialLawCapPress();

    /// Extract gas/oil capillary pressure value (Pg - Po) from current
    /// fluid state.
    double materialLawCapPressGasOil() const;

    /// Extract oil/water capillary pressure value (Po - Pw) from current
    /// fluid state.
    double materialLawCapPressOilWater() const;

    /// Extract gas/water capillary pressure value (Pg - Pw) from current
    /// fluid state.
    double materialLawCapPressGasWater() const;

    /// Predicate for whether specific phase has constant capillary pressure
    /// curve in current cell.
    ///
    /// \param[in] phaseIdx Phase.  Typically gas or water.
    ///
    /// \return Whether or not \p phaseIdx has constant capillary pressure
    /// curve in current cell.
    bool isConstCapPress(const PhaseIdx phaseIdx) const;

    /// Predicate for whether or not the G/O and O/W transition zones
    /// overlap in the current cell.
    ///
    /// This is the case when inverting the capillary pressure curves
    /// produces a negative oil saturation--i.e., when Sg + Sw > 1.
    bool isOverlappingTransition() const;

    /// Derive phase saturation value from simple depth consideration.
    ///
    /// Assumes that the pertinent capillary pressure curve is constant
    /// (typically zero) in the current cell--i.e., that there is a sharp
    /// interface between the two phases.
    ///
    /// \param[in] contactdepth Depth of relevant phase separation contact.
    ///
    /// \param[in] Position of phase in three-phase enumeration.  Typically
    ///    \code gasPos() \endcode or \code waterPos() \endcode.
    ///
    /// \param[in] isincr Whether the capillary pressure curve is normally
    ///    increasing as a function of phase saturation (e.g., Pcgo(Sg) = Pg
    ///    - Po) or if the curve is normally decreasing as a function of
    ///    increasing phase saturation (e.g., Pcow(Sw) = Po - Pw).  True for
    ///    capillary pressure functions that are normally increasing as a
    ///    function of phase saturation.
    ///
    /// \return Phase saturation.
    double fromDepthTable(const double   contactdepth,
                          const PhaseIdx phasePos,
                          const bool     isincr) const;

    /// Derive phase saturation by inverting non-constant capillary pressure
    /// curve.
    ///
    /// \param[in] pc Target capillary pressure value.
    ///
    /// \param[in] Position of phase in three-phase enumeration.  Typically
    ///    \code gasPos() \endcode or \code waterPos() \endcode.
    ///
    /// \param[in] isincr Whether the capillary pressure curve is normally
    ///    increasing as a function of phase saturation (e.g., Pcgo(Sg) = Pg
    ///    - Po) or if the curve is normally decreasing as a function of
    ///    increasing phase saturation (e.g., Pcow(Sw) = Po - Pw).  True for
    ///    capillary pressure functions that are normally increasing as a
    ///    function of phase saturation.
    ///
    /// \return Phase saturation at which capillary pressure attains target
    ///    value.
    double invertCapPress(const double   pc,
                          const PhaseIdx phasePos,
                          const bool     isincr) const;

    /// Position of oil in fluid system's three-phase enumeration.
    PhaseIdx oilPos() const
    {
        return FluidSystem::oilPhaseIdx;
    }

    /// Position of gas in fluid system's three-phase enumeration.
    PhaseIdx gasPos() const
    {
        return FluidSystem::gasPhaseIdx;
    }

    /// Position of water in fluid system's three-phase enumeration.
    PhaseIdx waterPos() const
    {
        return FluidSystem::waterPhaseIdx;
    }
};

// ===========================================================================

template <typename CellRange, typename Comm>
void verticalExtent(const CellRange&      cells,
                    const std::vector<std::pair<double, double>>& cellZMinMax,
                    const Comm& comm,
                    std::array<double,2>& span);

template <class Element>
std::pair<double,double> cellZMinMax(const Element& element);

} // namespace Details

namespace DeckDependent {

template<class FluidSystem,
         class Grid,
         class GridView,
         class ElementMapper,
         class CartesianIndexMapper>
class InitialStateComputer
{
    using Element = typename GridView::template Codim<0>::Entity;
public:
    template<class MaterialLawManager>
    InitialStateComputer(MaterialLawManager& materialLawManager,
                         const EclipseState& eclipseState,
                         const Grid& grid,
                         const GridView& gridView,
                         const CartesianIndexMapper& cartMapper,
                         const double grav,
                         const int num_pressure_points = 2000,
                         const bool applySwatInit = true);

    using Vec = std::vector<double>;
    using PVec = std::vector<Vec>; // One per phase.

    const Vec& temperature() const { return temperature_; }
    const Vec& saltConcentration() const { return saltConcentration_; }
    const Vec& saltSaturation() const { return saltSaturation_; }
    const PVec& press() const { return pp_; }
    const PVec& saturation() const { return sat_; }
    const Vec& rs() const { return rs_; }
    const Vec& rv() const { return rv_; }
    const Vec& rvw() const { return rvw_; }

private:
    template <class RMap>
    void updateInitialTemperature_(const EclipseState& eclState, const RMap& reg);

    template <class RMap>
    void updateInitialSaltConcentration_(const EclipseState& eclState, const RMap& reg);

    template <class RMap>
    void updateInitialSaltSaturation_(const EclipseState& eclState, const RMap& reg);

    void updateCellProps_(const GridView& gridView, 
                          const NumericalAquifers& aquifer);

    void applyNumericalAquifers_(const GridView& gridView,
                                 const NumericalAquifers& aquifer,
                                 const bool co2store_or_h2store);

    template<class RMap>
    void setRegionPvtIdx(const EclipseState& eclState, const RMap& reg);

    template <class RMap, class MaterialLawManager, class Comm>
    void calcPressSatRsRv(const RMap& reg,
                          const std::vector<EquilRecord>& rec,
                          MaterialLawManager& materialLawManager,
                          const Comm& comm,
                          const double grav);

    template <class CellRange, class EquilibrationMethod>
    void cellLoop(const CellRange&      cells,
                  EquilibrationMethod&& eqmethod);

    template <class CellRange, class PressTable, class PhaseSat>
    void equilibrateCellCentres(const CellRange&         cells,
                                const EquilReg&          eqreg,
                                const PressTable&        ptable,
                                PhaseSat&                psat);

    template <class CellRange, class PressTable, class PhaseSat>
    void equilibrateHorizontal(const CellRange&  cells,
                               const EquilReg&   eqreg,
                               const int         acc,
                               const PressTable& ptable,
                               PhaseSat&         psat);

    std::vector< std::shared_ptr<Miscibility::RsFunction> > rsFunc_;
    std::vector< std::shared_ptr<Miscibility::RsFunction> > rvFunc_;
    std::vector< std::shared_ptr<Miscibility::RsFunction> > rvwFunc_;
    using TabulatedFunction = Tabulated1DFunction<double>;
    std::vector<TabulatedFunction> tempVdTable_;
    std::vector<TabulatedFunction> saltVdTable_;
    std::vector<TabulatedFunction> saltpVdTable_;
    std::vector<int> regionPvtIdx_;
    Vec temperature_;
    Vec saltConcentration_;
    Vec saltSaturation_;
    PVec pp_;
    PVec sat_;
    Vec rs_;
    Vec rv_;
    Vec rvw_;
    const CartesianIndexMapper& cartesianIndexMapper_;
    Vec swatInit_;
    Vec cellCenterDepth_;
    std::vector<std::pair<double,double>> cellZSpan_;
    std::vector<std::pair<double,double>> cellZMinMax_;
    int num_pressure_points_;
};

} // namespace DeckDependent
} // namespace EQUIL
} // namespace Opm

#endif // OPM_INITSTATEEQUIL_HEADER_INCLUDED
