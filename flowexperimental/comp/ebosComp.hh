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
 * \brief The common settings for all ebos variants.
 */
#ifndef EBOSCOMP_HH
#define EBOSCOMP_HH


#include <opm/models/discretization/common/fvbaseproblem.hh>
#include <opm/models/utils/start.hh>

#include <opm/simulators/aquifers/BlackoilAquiferModel.hpp>

#include <opm/simulators/flow/FlowProblemComp.hpp>
#include <opm/simulators/flow/FlowProblemCompProperties.hpp>

#include <opm/simulators/linalg/ISTLSolver.hpp>
#include <opm/simulators/timestepping/EclTimeSteppingParams.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>

namespace Opm {
template <class TypeTag>
class EbosProblemComp;
}

namespace Opm {
template <class TypeTag>
class EbosProblemComp : public FlowProblemComp<TypeTag> //, public FvBaseProblem<TypeTag>
{
    typedef FlowProblemComp<TypeTag> ParentType; // FlowProblemComp
    using BaseType = GetPropType<TypeTag, Properties::BaseProblem>; // multiphase problem
public:
    void writeOutput(bool verbose = true)
    {
        OPM_TIMEBLOCK(problemWriteOutput);
        // use the generic code to prepare the output fields and to
        // write the desired VTK files.
        if (Parameters::get<TypeTag, Properties::EnableWriteAllSolutions>() || this->simulator().episodeWillBeOver()){
            BaseType::writeOutput(verbose);
        }
    }

    static void registerParameters()
    {
        // outputTypeTagInfo<ParentType>();
        // outputTypeTagInfo<BaseType>();
        ParentType::registerParameters();
        BaseType::registerParameters();

        BlackoilModelParameters<TypeTag>::registerParameters();
        Parameters::registerParam<TypeTag, Properties::EnableTerminalOutput>("Do *NOT* use!");
        Parameters::hideParam<TypeTag, Properties::DbhpMaxRel>();
        Parameters::hideParam<TypeTag, Properties::DwellFractionMax>();
        Parameters::hideParam<TypeTag, Properties::MaxResidualAllowed>();
        Parameters::hideParam<TypeTag, Properties::ToleranceMb>();
        Parameters::hideParam<TypeTag, Properties::ToleranceMbRelaxed>();
        Parameters::hideParam<TypeTag, Properties::ToleranceCnv>();
        Parameters::hideParam<TypeTag, Properties::ToleranceCnvRelaxed>();
        Parameters::hideParam<TypeTag, Properties::ToleranceWells>();
        Parameters::hideParam<TypeTag, Properties::ToleranceWellControl>();
        Parameters::hideParam<TypeTag, Properties::MaxWelleqIter>();
        Parameters::hideParam<TypeTag, Properties::UseMultisegmentWell>();
        Parameters::hideParam<TypeTag, Properties::TolerancePressureMsWells>();
        Parameters::hideParam<TypeTag, Properties::MaxPressureChangeMsWells>();
        Parameters::hideParam<TypeTag, Properties::MaxInnerIterMsWells>();
        Parameters::hideParam<TypeTag, Properties::MaxNewtonIterationsWithInnerWellIterations>();
        Parameters::hideParam<TypeTag, Properties::MaxInnerIterWells>();
        Parameters::hideParam<TypeTag, Properties::MaxSinglePrecisionDays>();
        Parameters::hideParam<TypeTag, Properties::MinStrictCnvIter>();
        Parameters::hideParam<TypeTag, Properties::MinStrictMbIter>();
        Parameters::hideParam<TypeTag, Properties::SolveWelleqInitially>();
        Parameters::hideParam<TypeTag, Properties::UpdateEquationsScaling>();
        Parameters::hideParam<TypeTag, Properties::UseUpdateStabilization>();
        Parameters::hideParam<TypeTag, Properties::MatrixAddWellContributions>();
        Parameters::hideParam<TypeTag, Properties::EnableTerminalOutput>();
    }

    // inherit the constructors
    using ParentType::FlowProblemComp;
};
}

#endif // EBOS_HH
