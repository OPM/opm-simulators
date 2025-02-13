/*
  Copyright 2024, SINTEF Digital

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

#ifndef OPM_COMP_WELL_HPP
#define OPM_COMP_WELL_HPP

#include <opm/models/utils/propertysystem.hh>

#include "CompWellEquations.hpp"
#include "CompWellInterface.hpp"
#include "CompWellPrimaryVariables.hpp"
#include "CompConnectionData.hpp"

#include <string>

namespace Opm
{

template <typename TypeTag>
class CompWell  : public CompWellInterface<TypeTag>
{
public:
    using Base = CompWellInterface<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using PrimaryVariables = CompWellPrimaryVariables<FluidSystem, Indices>;
    using WellEquations = CompWellEquations<Scalar, PrimaryVariables::numWellEq, Indices::numEq>;


    constexpr static unsigned num_comp = FluidSystem::numComponents;

    using EvalWell = typename PrimaryVariables::EvalWell;
    using BVectorWell = typename WellEquations::BVectorWell;

    // TODO: this can be a rate converter role later
    // currently, it has the surface densities for each phase and volume fractions for each phase
    // it is part of the secondary variables used in the assembling the well equations
    struct SurfaceConditons
    {
        static constexpr int num_phases = 2;
        std::array<EvalWell, num_phases> surface_densities_{0.};
        std::array<EvalWell, num_phases> volume_fractions_{0.};
        std::array<std::array<EvalWell, num_comp>, num_phases> mass_fractions_{0.};

        EvalWell density() const {
            EvalWell density = 0.;
            for (int i = 0; i < num_phases; ++i) {
                density += surface_densities_[i] * volume_fractions_[i];
            }
            return density;
        }

        // TODO: it looks like that we can have a concept of component mass density?
        EvalWell massFraction(int comp_idx) const {
            EvalWell mass = 0.;
            for (unsigned  p = 0; p < num_phases; ++p) {
                mass += surface_densities_[p] * volume_fractions_[p] * mass_fractions_[p][comp_idx];
            }
            return mass / density();
        }
    };

    CompWell(const Well& well,
             int index_of_well,
             const std::vector<CompConnectionData<Scalar>>& well_connection_data);

    void init() override;

    void calculateExplicitQuantities(const Simulator& simulator,
                                     const SingleCompWellState<Scalar>& well_state);

    void updatePrimaryVariables(const Simulator& simulator,
                                const SingleCompWellState<Scalar>& well_state) override;

    void updateSecondaryQuantities(const Simulator& simulator);

    void updatePrimaryVariableEvaluation(); // override;

    // TODO: control should be passed in later
    void assembleWellEq(const Simulator& simulator,
                        const double dt,
                        const Well::InjectionControls& inj_controls,
                        const Well::ProductionControls& prod_controls,
                        const SingleCompWellState<Scalar>& well_state);

    bool iterateWellEq(const Simulator& simulator,
                       const Scalar dt,
                       SingleCompWellState<Scalar>& well_state);

    void solveEqAndUpdateWellState(const Simulator& simulator,
                                   SingleCompWellState<Scalar>& well_state);


private:

    // primary variables
    PrimaryVariables primary_variables_;
    WellEquations well_equations_;

    // the following varialbes are temporary and remain to be cleaned up and re-organized
    // some are testing variables, and some are secondary variables might be kept
    // anyway, they are very rough prototype code for testing and will be changed
    const Scalar wellbore_volume_ {21.6}; // m^3, it is rather big, will come with different design when the working flow is established

    std::array<EvalWell, num_comp> mass_fractions_{0.};
    EvalWell fluid_density_{0.};
    // the original mass for each component in wellbore
    std::array<Scalar, num_comp> component_masses_{0.};
    // the new mass for each component in wellbore, derived from the primary variables
    std::array<EvalWell, num_comp> new_component_masses_{0.};
    // quantities used to calculate the quantities under the surface conditions
    SurfaceConditons surface_conditions_;

    // following are some secondary property or variables to be used for later
    void calculateSingleConnectionRate(const Simulator& simulator,
                                       std::vector<EvalWell>& cq_s) const;

    void updateTotalMass();
    // TODO: a better name
    void updateSurfaceQuantities(const Simulator& simulator);

    void getMoblity(const Simulator& simulator,
                    int connection_idx,
                    std::vector<EvalWell>& mob) const;


    // TODO: the following assembling functions will be moved to a separate assmeble class
    void assembleSourceTerm(const Scalar dt);

    void updatePrimaryVariablesNewton(const BVectorWell& dwells);
};

} // end of namespace Opm

#include "CompWell_impl.hpp"

#endif // OPM_COMPOSITIONAL_WELL_MODEL_HPP