/*
  Copyright 2020, NORCE AS

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
#include <config.h>

#include <ebos/eclpolyhedralgridvanguard.hh>

#include <opm/grid/polyhedralgrid.hh>
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
#include <opm/simulators/flow/Main.hpp>

// these are not explicitly instanced in library
#include <ebos/collecttoiorank_impl.hh>
#include <ebos/eclgenericproblem_impl.hh>
#include <ebos/eclgenericthresholdpressure_impl.hh>
#include <ebos/eclgenerictracermodel_impl.hh>
#include <ebos/ecltransmissibility_impl.hh>
#include <ebos/eclgenericwriter_impl.hh>
#include <ebos/equil/initstateequil_impl.hh>
#include <opm/simulators/utils/GridDataOutput_impl.hpp>

namespace Opm {
namespace Properties {
    namespace TTag {
        struct FlowProblemPoly {
            using InheritsFrom = std::tuple<FlowProblem>;
        };
    }

    template<class TypeTag>
    struct Linearizer<TypeTag, TTag::FlowProblemPoly> { using type = TpfaLinearizer<TypeTag>; };

    template<class TypeTag>
    struct LocalResidual<TypeTag, TTag::FlowProblemPoly> { using type = BlackOilLocalResidualTPFA<TypeTag>; };

    template<class TypeTag>
    struct EnableDiffusion<TypeTag, TTag::FlowProblemPoly> { static constexpr bool value = false; };

    template<class TypeTag>
    struct Grid<TypeTag, TTag::FlowProblemPoly> {
        using type = Dune::PolyhedralGrid<3, 3>;
    };
    template<class TypeTag>
    struct EquilGrid<TypeTag, TTag::FlowProblemPoly> {
        //using type = Dune::CpGrid;
        using type = GetPropType<TypeTag, Properties::Grid>;
    };

    template<class TypeTag>
    struct Vanguard<TypeTag, TTag::FlowProblemPoly> {
        using type = Opm::EclPolyhedralGridVanguard<TypeTag>;
    };
}
}

int main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::FlowProblemPoly;
    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    auto ret = mainObject->runStatic<TypeTag>();
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}
