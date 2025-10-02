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
// By default, dune-ALUGrid uses Space-Filling Curve (SFC) ordering
// for cells to optimize data access patterns for adaptive mesh refinements/coarsening.
// However, if you want to use Cartesian ordering, as is used in OPM, you
// can switch to it by defining the macro #define DISABLE_SFC_ORDERING 1.
// This will change the default cell order to Cartesian.
// Note that this option is not available for pre-built or installed versions of dune-ALUGrid.
// To enable changig to Cartesian ordering, you will need to rebuild dune-ALUGrid from source, ensuring
// the build configuration allows disabling SFC ordering from OPM.
// For more details, refer to the files gridfactory.hh and aluinline.hh located in the dune-alugrid/3d/

#include <dune/alugrid/grid.hh>
#include <opm/simulators/flow/Main.hpp>

// for equilgrid in writer
// need to include this before eclgenericwriter_impl.hh due to specializations.
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/GridHelpers.hpp>

// these are not explicitly instanced in library
#include <opm/simulators/flow/AluGridVanguard.hpp>
#include <opm/simulators/flow/CollectDataOnIORank_impl.hpp>
#include <opm/simulators/flow/EclGenericWriter_impl.hpp>
#include <opm/simulators/flow/FlowGenericProblem_impl.hpp>
#include <opm/simulators/flow/GenericThresholdPressure_impl.hpp>
#include <opm/simulators/flow/GenericTracerModel_impl.hpp>
#include <opm/simulators/flow/Transmissibility_impl.hpp>
#include <opm/simulators/flow/equil/InitStateEquil_impl.hpp>
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
    using type = Opm::AluGridVanguard<TypeTag>;
};
}

template<>
class SupportsFaceTag<Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>>
    : public std::bool_constant<true>
{};

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
