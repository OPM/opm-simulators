/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2016 - 2017 IRIS AS.

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


#ifndef OPM_STANDARDWELL_CONNECTIONS_HEADER_INCLUDED
#define OPM_STANDARDWELL_CONNECTIONS_HEADER_INCLUDED

#include <opm/simulators/wells/StandardWellPrimaryVariables.hpp>

#include <array>
#include <functional>
#include <tuple>
#include <variant>
#include <vector>

namespace Opm
{

class DeferredLogger;
enum class Phase;
template<typename FluidSystem, typename Indices> class WellInterfaceIndices;
template<typename Scalar, typename IndexTraits> class WellState;
template<class Scalar> class PerfData;

template<class FluidSystem, class Indices>
class StandardWellConnections
{
public:
    using Scalar = typename FluidSystem::Scalar;
    using IndexTraits = typename FluidSystem::IndexTraitsType;
    explicit StandardWellConnections(const WellInterfaceIndices<FluidSystem, Indices>& well);

    struct Properties
    {
        std::vector<Scalar> b_perf{};
        std::vector<Scalar> rsmax_perf{};
        std::vector<Scalar> rvmax_perf{};
        std::vector<Scalar> rvwmax_perf{};
        std::vector<Scalar> rswmax_perf{};
        std::vector<Scalar> surf_dens_perf{};
    };

    struct PressurePropertyFunctions
    {
        std::function<Scalar(int,int)> getTemperature{};
        std::function<Scalar(int)>     getSaltConcentration{};
        std::function<int(int)>        pvtRegionIdx{};
        std::function<Scalar(int)>     solventInverseFormationVolumeFactor{};
        std::function<Scalar(int)>     solventRefDensity{};
    };

    struct DensityPropertyFunctions
    {
        std::function<void(int, const std::vector<int>&, std::vector<Scalar>&)> mobility{};
        std::function<void(int, const std::vector<int>&, std::vector<Scalar>&)> densityInCell{};
    };

    Properties
    computePropertiesForPressures(const WellState<Scalar, IndexTraits>&         well_state,
                                  const PressurePropertyFunctions& propFunc) const;

    //! \brief Compute connection properties (densities, pressure drop, ...)
    void computeProperties(const bool                      stop_or_zero_rate_target,
                           const WellState<Scalar, IndexTraits>&        well_state,
                           const DensityPropertyFunctions& prop_func,
                           const Properties&               props,
                           DeferredLogger&                 deferred_logger);

    //! \brief Returns density for specific perforation/connection.
    //!
    //! \param[in] i Connection index
    //!
    //! \return Mixture density at connection \p i.
    Scalar rho(const typename std::vector<Scalar>::size_type i) const
    {
        return (i < this->perf_densities_.size())
            ? this->perf_densities_[i]
            : 0.0;
    }

    //! \brief Returns pressure drop for a given perforation.
    Scalar pressure_diff(const unsigned perf) const
    { return perf_pressure_diffs_[perf]; }

    using Eval = typename WellInterfaceIndices<FluidSystem, Indices>::Eval;
    using EvalWell = typename StandardWellPrimaryVariables<FluidSystem, Indices>::EvalWell;

    Eval connectionRateBrine(Scalar& rate,
                             const Scalar vap_wat_rate,
                             const std::vector<EvalWell>& cq_s,
                             const std::variant<Scalar,EvalWell>& saltConcentration) const;

    Eval connectionRateFoam(const std::vector<EvalWell>& cq_s,
                            const std::variant<Scalar,EvalWell>& foamConcentration,
                            const Phase transportPhase,
                            DeferredLogger& deferred_logger) const;

    std::tuple<Eval,EvalWell>
    connectionRatePolymer(Scalar& rate,
                          const std::vector<EvalWell>& cq_s,
                          const std::variant<Scalar,EvalWell>& polymerConcentration) const;

    Eval connectionRateBioeffects(Scalar& rate,
                                  const Scalar vap_wat_rate,
                                  const std::vector<EvalWell>& cq_s,
                                  const std::variant<Scalar,EvalWell>& microbialConcentration) const;

    std::tuple<Eval,Eval,Eval>
    connectionRatesMICP(Scalar& rate_m,
                        Scalar& rate_o,
                        Scalar& rate_u,
                        const std::vector<EvalWell>& cq_s,
                        const std::variant<Scalar,EvalWell>& microbialConcentration,
                        const std::variant<Scalar,EvalWell>& oxygenConcentration,
                        const std::variant<Scalar,EvalWell>& ureaConcentration) const;

    std::tuple<Eval,EvalWell>
    connectionRatezFraction(Scalar& rate,
                            const Scalar dis_gas_rate,
                            const std::vector<EvalWell>& cq_s,
                            const std::variant<Scalar, std::array<EvalWell,2>>& solventConcentration) const;

private:
    void computePressureDelta();

    // TODO: not total sure whether it is a good idea to put this function here
    // the major reason to put here is to avoid the usage of Wells struct
    void computeDensities(const std::vector<Scalar>& perfComponentRates,
                          const Properties& props,
                          DeferredLogger& deferred_logger);

    void computeDensitiesForStoppedProducer(const DensityPropertyFunctions& prop_func);

    std::vector<Scalar>
    calculatePerforationOutflow(const std::vector<Scalar>& perfComponentRates) const;

    void initialiseConnectionMixture(const int                  num_comp,
                                     const int                  perf,
                                     const std::vector<Scalar>& q_out_perf,
                                     const std::vector<Scalar>& currentMixture,
                                     std::vector<Scalar>&       previousMixture) const;

    std::vector<Scalar>
    copyInPerforationRates(const Properties&       props,
                           const PerfData<Scalar>& perf_data) const;

    const WellInterfaceIndices<FluidSystem, Indices>& well_; //!< Reference to well interface

    std::vector<Scalar> perf_densities_; //!< densities of the fluid in each perforation
    std::vector<Scalar> perf_pressure_diffs_; //!< // pressure drop between different perforations
};

}

#endif // OPM_STANDARDWELL_CONNECTIONS_HEADER_INCLUDED
