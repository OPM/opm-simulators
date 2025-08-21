/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 - 2017 Statoil ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 IRIS AS

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
#ifndef OPM_BLACKOILWELLMODEL_NLDD_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_NLDD_HEADER_INCLUDED

#include <opm/grid/utility/SparseTable.hpp>

#include <opm/simulators/timestepping/ConvergenceReport.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>

#include <opm/simulators/wells/BlackoilWellModel.hpp>

#include <map>
#include <string>
#include <vector>

namespace Opm {

template<typename Scalar, typename IndexTraits>
class BlackoilWellModelNlddGeneric
{
public:
    std::vector<Scalar> getPrimaryVarsDomain(const int domainIdx) const;
    void setPrimaryVarsDomain(const int domainIdx, const std::vector<Scalar>& vars);

    const SparseTable<int>& well_local_cells() const
    { return well_local_cells_; }

    const std::map<std::string, int>& well_domain() const
    { return well_domain_; }

protected:
    BlackoilWellModelNlddGeneric(BlackoilWellModelGeneric<Scalar, IndexTraits>& model)
        : genWellModel_(model)
    {}

    void calcDomains(const std::vector<const SubDomainIndices*>& domains);

private:
    void logDomains() const;

    void findWellDomains(const std::vector<const SubDomainIndices*>& domains);

    void calcLocalIndices(const std::vector<const SubDomainIndices*>& domains);

    BlackoilWellModelGeneric<Scalar, IndexTraits>& genWellModel_;

    // Keep track of the domain of each well
    std::map<std::string, int> well_domain_{};

    // Store the local index of the wells perforated cells in the domain
    SparseTable<int> well_local_cells_{};
};

/// Class for handling the blackoil well model in a NLDD solver.
template<typename TypeTag>
class BlackoilWellModelNldd :
    public BlackoilWellModelNlddGeneric<GetPropType<TypeTag, Properties::Scalar>,
                                        typename GetPropType<TypeTag, Properties::FluidSystem>::IndexTraitsType>
{
public:
    // ---------      Types      ---------
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PressureMatrix = typename BlackoilWellModel<TypeTag>::PressureMatrix;
    using BVector = typename BlackoilWellModel<TypeTag>::BVector;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using IndexTraits = typename FluidSystem::IndexTraitsType;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    using Domain = SubDomain<Grid>;

    BlackoilWellModelNldd(BlackoilWellModel<TypeTag>& model)
        : BlackoilWellModelNlddGeneric<Scalar, IndexTraits>(model)
        , wellModel_(model)
    {}

    void addWellPressureEquations(PressureMatrix& jacobian,
                                  const BVector& weights,
                                  const bool use_well_weights,
                                  const int domainIndex) const;

    // prototype for assemble function for ASPIN solveLocal()
    // will try to merge back to assemble() when done prototyping
    void assemble(const int iterationIdx,
                  const double dt,
                  const Domain& domain);

    void updateWellControls(DeferredLogger& deferred_logger,
                            const Domain& domain);

    void setupDomains(const std::vector<Domain>& domains);

    // Check if well equations are converged locally.
    ConvergenceReport getWellConvergence(const Domain& domain,
                                         const std::vector<Scalar>& B_avg,
                                         DeferredLogger& local_deferredLogger) const;

    // using the solution x to recover the solution xw for wells and applying
    // xw to update Well State
    void recoverWellSolutionAndUpdateWellState(const BVector& x,
                                               const int domainIdx);

    // Get number of wells on this rank
    int numLocalWells() const 
    {
        return wellModel_.numLocalWells(); 
    }

    // Get number of wells on this rank
    int numLocalWellsEnd() const 
    {
        return wellModel_.numLocalWellsEnd(); 
    }

private:
    BlackoilWellModel<TypeTag>& wellModel_;

    void assembleWellEq(const double dt,
                        const Domain& domain,
                        DeferredLogger& deferred_logger);

    // These members are used to avoid reallocation in specific functions
    // instead of using local variables.
    // Their state is not relevant between function calls, so they can
    // (and must) be mutable, as the functions using them are const.
    mutable BVector x_local_;

};

} // namespace Opm

#include "BlackoilWellModelNldd_impl.hpp"

#endif
