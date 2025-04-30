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
 * \copydoc Opm::FvBaseDiscretization
 */
#ifndef EWOMS_FV_BASE_DISCRETIZATION_FEMADAPT_HH
#define EWOMS_FV_BASE_DISCRETIZATION_FEMADAPT_HH

#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/misc/capabilities.hh>
#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/space/common/restrictprolongtuple.hh>

#include <opm/models/discretization/common/fvbasediscretization.hh>

#include <memory>
#include <stdexcept>

namespace Opm {

template<class TypeTag>
class FvBaseDiscretizationFemAdapt;

namespace Properties {

template<class TypeTag>
struct BaseDiscretizationType<TypeTag, TTag::FvBaseDiscretization>
{ using type = FvBaseDiscretizationFemAdapt<TypeTag>; };

template<class TypeTag>
struct DiscreteFunction<TypeTag, TTag::FvBaseDiscretization>
{
    using DiscreteFunctionSpace  = GetPropType<TypeTag, Properties::DiscreteFunctionSpace>;
    using PrimaryVariables  = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using type = Dune::Fem::ISTLBlockVectorDiscreteFunction<DiscreteFunctionSpace, PrimaryVariables>;
};

} // namespace Properties

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief The base class for the finite volume discretization schemes.
 */

template <class TypeTag>
class FvBaseDiscretizationFemAdapt : public FvBaseDiscretization<TypeTag>
{
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using ParentType = FvBaseDiscretization<TypeTag>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

    static constexpr unsigned historySize = getPropValue<TypeTag, Properties::TimeDiscHistorySize>();

    using DiscreteFunctionSpace = GetPropType<TypeTag, Properties::DiscreteFunctionSpace>;
    // discrete function storing solution data
    using DiscreteFunction = Dune::Fem::ISTLBlockVectorDiscreteFunction<DiscreteFunctionSpace, PrimaryVariables>;

    // problem restriction and prolongation operator for adaptation
    using ProblemRestrictProlongOperator = typename Problem::RestrictProlongOperator;

    // discrete function restriction and prolongation operator for adaptation
    using DiscreteFunctionRestrictProlong = Dune::Fem::RestrictProlongDefault<DiscreteFunction>;
    using RestrictProlong
        = Dune::Fem::RestrictProlongTuple<DiscreteFunctionRestrictProlong, ProblemRestrictProlongOperator>;

    // adaptation classes
    using AdaptationManager = Dune::Fem::AdaptationManager<Grid, RestrictProlong>;

public:
    template<class Serializer>
    struct SerializeHelper
    {
        template<class SolutionType>
        static void serializeOp(Serializer& serializer,
                                SolutionType& solution)
        {
            for (auto& sol : solution) {
                serializer(sol->blockVector());
            }
        }
    };

    explicit FvBaseDiscretizationFemAdapt(Simulator& simulator)
        : ParentType(simulator)
        , space_(simulator.vanguard().gridPart())
    {
        if (this->enableGridAdaptation_ && !Dune::Fem::Capabilities::isLocallyAdaptive<Grid>::v) {
            throw std::invalid_argument("Grid adaptation enabled, but chosen Grid is not capable"
                                        " of adaptivity");
        }

        for (unsigned timeIdx = 0; timeIdx < historySize; ++timeIdx) {
            this->solution_[timeIdx] = std::make_unique<DiscreteFunction>("solution", space_);
        }
    }

    void adaptGrid()
    {
        // adapt the grid if enabled and if all dependencies are available
        // adaptation is only done if markForGridAdaptation returns true
        if (this->enableGridAdaptation_) {
            // check if problem allows for adaptation and cells were marked
            if (this->simulator_.problem().markForGridAdaptation()) {
                // adapt the grid and load balance if necessary
                adaptationManager().adapt();

                // if the grid has potentially changed, we need to re-create the
                // supporting data structures.
                this->elementMapper_.update(this->gridView_);
                this->vertexMapper_.update(this->gridView_);
                this->resetLinearizer();
                // this is a bit hacky because it supposes that Problem::finishInit()
                // works fine multiple times in a row.
                //
                // TODO: move this to Problem::gridChanged()
                this->finishInit();

                // notify the problem that the grid has changed
                //
                // TODO: come up with a mechanism to access the unadapted data structures
                // outside of the problem (i.e., grid, mappers, solutions)
                this->simulator_.problem().gridChanged();

                // notify the modules for visualization output
                for (auto& module : this->outputModules_) {
                    module->allocBuffers();
                }
            }
        }
    }

    AdaptationManager& adaptationManager()
    {
        if (!adaptationManager_) {
            // create adaptation objects here, because when doing so in constructor
            // problem is not yet intialized, aka seg fault
            restrictProlong_ =
                std::make_unique<RestrictProlong>(DiscreteFunctionRestrictProlong(*(this->solution_[/*timeIdx=*/0])),
                                                  this->simulator_.problem().restrictProlongOperator());
            adaptationManager_ =
                std::make_unique<AdaptationManager>(this->simulator_.vanguard().grid(), *restrictProlong_);
        }
        return *adaptationManager_;
    }

private:
    DiscreteFunctionSpace space_;
    std::unique_ptr<RestrictProlong> restrictProlong_;
    std::unique_ptr<AdaptationManager> adaptationManager_;
};

} // namespace Opm

#endif
