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
#include "config.h"
#include <opm/simulators/flow/Main.hpp>
#include <dune/alugrid/grid.hh>
#include <ebos/eclalugridvanguard.hh>
namespace Opm {
namespace Properties {
namespace TTag {
struct EclFlowProblemAlugrid {
    using InheritsFrom = std::tuple<EclFlowProblem>;
};
}

template<class TypeTag>
struct Grid<TypeTag, TTag::EclFlowProblemAlugrid> {
    static const int dim = 3;
    using type = Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming,Dune::ALUGridMPIComm>;
};
// alugrid need cp grid as equilgrid
template<class TypeTag>
struct EquilGrid<TypeTag, TTag::EclFlowProblemAlugrid> {
    using type = Dune::CpGrid;
};
template<class TypeTag>
struct Vanguard<TypeTag, TTag::EclFlowProblemAlugrid> {
    using type = Opm::EclAluGridVanguard<TypeTag>;
};
template<class TypeTag>
struct EclEnableAquifers<TypeTag, TTag::EclFlowProblemAlugrid> {
    static constexpr bool value = false;
};
}
}
int main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::EclFlowProblemAlugrid;
    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    auto ret = mainObject->runStatic<TypeTag>();
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}
