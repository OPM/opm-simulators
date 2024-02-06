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

#include <dune/alugrid/grid.hh>
#include <ebos/eclalugridvanguard.hh>
#include <opm/simulators/flow/Main.hpp>

// for equilgrid in writer
// need to include this before eclgenericwriter_impl.hh due to specializations.
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/GridHelpers.hpp>

// these are not explicitly instanced in library
#include <ebos/collecttoiorank_impl.hh>
#include <ebos/eclgenericproblem_impl.hh>
#include <ebos/eclgenericthresholdpressure_impl.hh>
#include <ebos/eclgenerictracermodel_impl.hh>
#include <ebos/eclgenericwriter_impl.hh>
#include <ebos/ecltransmissibility_impl.hh>
#include <ebos/equil/initstateequil_impl.hh>
#include <opm/simulators/utils/GridDataOutput_impl.hpp>

namespace Opm {
namespace Properties {
namespace TTag {
struct FlowProblemAlugrid {
    using InheritsFrom = std::tuple<FlowProblem>;
};
}

template<class TypeTag>
struct Grid<TypeTag, TTag::FlowProblemAlugrid> {
    static const int dim = 3;
#if HAVE_MPI
     using type = Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming,Dune::ALUGridMPIComm>;
#else
     using type = Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming, Dune::ALUGridNoComm>;
#endif
};

// alugrid need cp grid as equilgrid
template<class TypeTag>
struct EquilGrid<TypeTag, TTag::FlowProblemAlugrid> {
    using type = Dune::CpGrid;
};
template<class TypeTag>
struct Vanguard<TypeTag, TTag::FlowProblemAlugrid> {
    using type = Opm::EclAluGridVanguard<TypeTag>;
};
}
}

int main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::FlowProblemAlugrid;
    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    auto ret = mainObject->runStatic<TypeTag>();
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}
