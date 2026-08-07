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
 * \brief Black-oil specific routines for hydrostatic equilibrium-based
 *        initialisation, including dissolved gas (Rs), vaporized oil (Rv),
 *        and vaporized water (Rvw).
 */
#ifndef OPM_INIT_STATE_EQUIL_HPP
#define OPM_INIT_STATE_EQUIL_HPP

#include <opm/simulators/flow/equil/InitStateEquilBase.hpp>

namespace Opm {

namespace EQUIL {

template<class Scalar> class EquilReg;
namespace Miscibility { template<class Scalar> class RsFunction; }

namespace DeckDependent {

template<class FluidSystem,
         class Grid,
         class GridView,
         class ElementMapper,
         class CartesianIndexMapper>
class InitialStateComputer
    : public InitialStateComputerBase<FluidSystem, Grid, GridView,
                                      ElementMapper, CartesianIndexMapper>
{
    using Base = InitialStateComputerBase<FluidSystem, Grid, GridView,
                                          ElementMapper, CartesianIndexMapper>;
    using Scalar = typename FluidSystem::Scalar;

public:
    using typename Base::Vec;
    using typename Base::PVec;

    template<class MaterialLawManager>
    InitialStateComputer(MaterialLawManager& materialLawManager,
                         const EclipseState& eclipseState,
                         const Grid& grid,
                         const GridView& gridView,
                         const CartesianIndexMapper& cartMapper,
                         const Scalar grav,
                         const int num_pressure_points = 2000,
                         const bool applySwatInit = true);

    using Base::temperature;
    using Base::saltConcentration;
    using Base::saltSaturation;
    using Base::press;
    using Base::saturation;

    const Vec& rs() const { return rs_; }
    const Vec& rv() const { return rv_; }
    const Vec& rvw() const { return rvw_; }

private:
    template <class RMap, class MaterialLawManager, class Comm>
    void calcPressSatRsRv(const RMap& reg,
                          const std::vector<EquilRecord>& rec,
                          MaterialLawManager& materialLawManager,
                          const GridView& gridView,
                          const Comm& comm,
                          const Scalar grav);

    template <class CellRange, class EquilibrationMethod>
    void cellLoop(const CellRange&      cells,
                  EquilibrationMethod&& eqmethod);

    template <class CellRange, class PressTable, class PhaseSat>
    void equilibrateCellCentres(const CellRange&        cells,
                                const EquilReg<Scalar>& eqreg,
                                const PressTable&       ptable,
                                PhaseSat&               psat);

    template <class CellRange, class PressTable, class PhaseSat>
    void equilibrateHorizontal(const CellRange&        cells,
                               const EquilReg<Scalar>& eqreg,
                               const int               acc,
                               const PressTable&       ptable,
                               PhaseSat&               psat);

     template<class CellRange, class PressTable, class PhaseSat>
     void equilibrateTiltedFaultBlock(const CellRange& cells,
                            const EquilReg<Scalar>& eqreg,
                            const GridView& gridView, const int numLevels,
                            const PressTable& ptable, PhaseSat& psat);

     template<class CellRange, class PressTable, class PhaseSat>
     void equilibrateTiltedFaultBlockSimple(const CellRange& cells,
                           const EquilReg<Scalar>& eqreg,
                           const GridView& gridView, const int numLevels,
                           const PressTable& ptable, PhaseSat& psat);

    std::vector< std::shared_ptr<Miscibility::RsFunction<Scalar>> > rsFunc_;
    std::vector< std::shared_ptr<Miscibility::RsFunction<Scalar>> > rvFunc_;
    std::vector< std::shared_ptr<Miscibility::RsFunction<Scalar>> > rvwFunc_;
    Vec rs_;
    Vec rv_;
    Vec rvw_;
};

} // namespace DeckDependent
} // namespace EQUIL
} // namespace Opm

#endif // OPM_INIT_STATE_EQUIL_HPP
