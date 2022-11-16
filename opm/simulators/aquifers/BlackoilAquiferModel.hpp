/*
  File adapted from BlackoilWellModel.hpp

  Copyright 2017 TNO - Heat Transfer & Fluid Dynamics, Modelling & Optimization of the Subsurface
  Copyright 2017 Statoil ASA.

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
*/


#ifndef OPM_BLACKOILAQUIFERMODEL_HEADER_INCLUDED
#define OPM_BLACKOILAQUIFERMODEL_HEADER_INCLUDED

#include <ebos/eclbaseaquifermodel.hh>

#include <opm/input/eclipse/EclipseState/Aquifer/Aquancon.hpp>
#include <opm/input/eclipse/EclipseState/Aquifer/AquiferCT.hpp>
#include <opm/input/eclipse/EclipseState/Aquifer/Aquifetp.hpp>

#include <opm/output/data/Aquifer.hpp>

#include <opm/simulators/aquifers/AquiferCarterTracy.hpp>
#include <opm/simulators/aquifers/AquiferFetkovich.hpp>
#include <opm/simulators/aquifers/AquiferNumerical.hpp>

#include <opm/grid/CpGrid.hpp>
#ifdef USE_POLYHEDRALGRID
#include <opm/grid/polyhedralgrid.hh>
#endif
#if USE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <opm/material/densead/Math.hpp>

#include <vector>
#include <type_traits>

namespace Opm
{

template<class Grid>
class SupportsFaceTag
    : public std::bool_constant<false>
{};


template<>
class SupportsFaceTag<Dune::CpGrid>
    : public std::bool_constant<true>
{};


#ifdef USE_POLYHEDRALGRID
template<>
class SupportsFaceTag<Dune::PolyhedralGrid<3, 3>>
    : public std::bool_constant<true>
{};
#endif

#if USE_ALUGRID
template<>
class SupportsFaceTag<Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>>
    : public std::bool_constant<true>
{};
#endif


/// Class for handling the blackoil well model.
template <typename TypeTag>
class BlackoilAquiferModel
{
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;


public:
    explicit BlackoilAquiferModel(Simulator& simulator);

    void initialSolutionApplied();
    void initFromRestart(const data::Aquifers& aquiferSoln);

    void beginEpisode();
    void beginTimeStep();
    void beginIteration();
    // add the water rate due to aquifers to the source term.
    template <class Context>
    void addToSource(RateVector& rates, const Context& context, unsigned spaceIdx, unsigned timeIdx) const;
    void addToSource(RateVector& rates, unsigned globalSpaceIdx, unsigned timeIdx) const;
    void endIteration();
    void endTimeStep();
    void endEpisode();

    data::Aquifers aquiferData() const;

    template <class Restarter>
    void serialize(Restarter& res);

    template <class Restarter>
    void deserialize(Restarter& res);

protected:
    // ---------      Types      ---------
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    Simulator& simulator_;

    std::vector<std::unique_ptr<AquiferInterface<TypeTag>>> aquifers;

    // This initialization function is used to connect the parser objects with the ones needed by AquiferCarterTracy
    void init();
};


} // namespace Opm

#include "BlackoilAquiferModel_impl.hpp"

#endif
