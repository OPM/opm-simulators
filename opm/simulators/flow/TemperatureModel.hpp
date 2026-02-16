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
/**
 * \file
 *
 * \copydoc Opm::TemperatureModel
 */
#ifndef OPM_TEMPERATURE_MODEL_HPP
#define OPM_TEMPERATURE_MODEL_HPP

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/models/parallel/gridcommhandles.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/BlackoilModelParameters.hpp>
#include <opm/simulators/flow/GenericTemperatureModel.hpp>
#include <opm/simulators/linalg/findOverlapRowsAndColumns.hpp>
#include <opm/simulators/aquifers/AquiferGridUtils.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace Opm::Properties {

template<class TypeTag, class MyTypeTag>
struct EnableTemperatureModel {
    using type = UndefinedProperty;
};

} // namespace Opm::Properties

namespace Opm {

template<typename Scalar, typename IndexTraits> class WellState;

/*!
 * \ingroup BlackOilSimulator
 *
 * \brief A class which handles sequential implicit solution of the energy equation as specified in by TEMP
 */
template <class TypeTag, bool enableTempV = getPropValue<TypeTag, Properties::EnergyModuleType>() == EnergyModules::SequentialImplicitThermal >
class TemperatureModel : public GenericTemperatureModel<GetPropType<TypeTag, Properties::Grid>,
                                              GetPropType<TypeTag, Properties::GridView>,
                                              GetPropType<TypeTag, Properties::DofMapper>,
                                              GetPropType<TypeTag, Properties::Stencil>,
                                              GetPropType<TypeTag, Properties::FluidSystem>,
                                              GetPropType<TypeTag, Properties::Scalar>>
{
    using BaseType = GenericTemperatureModel<GetPropType<TypeTag, Properties::Grid>,
                                        GetPropType<TypeTag, Properties::GridView>,
                                        GetPropType<TypeTag, Properties::DofMapper>,
                                        GetPropType<TypeTag, Properties::Stencil>,
                                        GetPropType<TypeTag, Properties::FluidSystem>,
                                        GetPropType<TypeTag, Properties::Scalar>>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Stencil = GetPropType<TypeTag, Properties::Stencil>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using TemperatureEvaluation = DenseAd::Evaluation<Scalar,1>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using EnergyModule = BlackOilEnergyModule<TypeTag>;
    using IndexTraits = typename FluidSystem::IndexTraitsType;
    using WellStateType = WellState<Scalar, IndexTraits>;

    using EnergyMatrix = typename BaseType::EnergyMatrix;
    using EnergyVector = typename BaseType::EnergyVector;

    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { numPhases = FluidSystem::numPhases };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };

public:
    explicit TemperatureModel(Simulator& simulator)
        : BaseType(simulator.vanguard().gridView(),
                   simulator.vanguard().eclState(),
                   simulator.vanguard().cartesianIndexMapper(),
                   simulator.model().dofMapper())
        , simulator_(simulator)
    {}

    void init()
    {
        const unsigned int numCells = simulator_.model().numTotalDof();
        this->doInit(numCells);

        if (!this->doTemp())
            return;

        // we need the storage term at start of the iteration (timeIdx = 1)
        storage1_.resize(numCells);

        // set the initial temperature
        for (unsigned globI = 0; globI < numCells; ++globI) {
            this->temperature_[globI] = simulator_.problem().initialFluidState(globI).temperature(0);
        }
        // keep a copy of the intensive quantities to simplify the update during
        // the newton iterations
        intQuants_.resize(numCells);

        // find and store the overlap cells
        const auto& elemMapper = simulator_.model().elementMapper();
        detail::findOverlapAndInterior(simulator_.vanguard().grid(), elemMapper, overlapRows_, interiorRows_);
    }

    void beginTimeStep()
    {
        if (!this->doTemp()) {
            return;
        }

        // We copy the intensive quantities here to make it possible to update them
        const unsigned int numCells = simulator_.model().numTotalDof();
        for (unsigned globI = 0; globI < numCells; ++globI) {
            intQuants_[globI] = simulator_.model().intensiveQuantities(globI, /*timeIdx*/ 0);
            intQuants_[globI].updateTemperature_(simulator_.problem(), globI, /*timeIdx*/ 0);
            intQuants_[globI].updateEnergyQuantities_(simulator_.problem(), globI, /*timeIdx*/ 0);
        }
        updateStorageCache();

        const int nw = simulator_.problem().wellModel().wellState().numWells();
        this->energy_rates_.resize(nw, 0.0);
    }

    /*!
     * \brief Informs the temperature model that a time step has just been finished.
     */
    void endTimeStep(WellStateType& wellState)
    {
        if (!this->doTemp()) {
            return;
        }

        // We copy the intensive quantities here to make it possible to update them
        const unsigned int numCells = simulator_.model().numTotalDof();
        for (unsigned globI = 0; globI < numCells; ++globI) {
            intQuants_[globI] = simulator_.model().intensiveQuantities(globI, /*timeIdx*/ 0);
            intQuants_[globI].updateTemperature_(simulator_.problem(), globI, /*timeIdx*/ 0);
            intQuants_[globI].updateEnergyQuantities_(simulator_.problem(), globI, /*timeIdx*/ 0);
        }
        advanceTemperatureFields();

        // update energy_rates
        const int nw = wellState.numWells();
        for (auto wellID = 0*nw; wellID < nw; ++wellID) {
            auto& ws = wellState.well(wellID);
            ws.energy_rate = this->energy_rates_[wellID];
        }
    }

    /*!
     * \brief This method writes the complete state of all temperature
     *        to the hard disk.
     */
    template <class Restarter>
    void serialize(Restarter&)
    { /* not implemented */ }

    /*!
     * \brief This method restores the complete state of the temperature
     *        from disk.
     *
     * It is the inverse of the serialize() method.
     */
    template <class Restarter>
    void deserialize(Restarter&)
    { /* not implemented */ }

protected:
    void updateStorageCache()
    {
        // we need the storage term at start of the iteration (timeIdx = 1)
        const unsigned int numCells = simulator_.model().numTotalDof();
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (unsigned globI = 0; globI < numCells; ++globI) {
            Scalar storage = 0.0;
            computeStorageTerm(globI, storage);
            storage1_[globI] = storage;
        }
    }

    void advanceTemperatureFields()
    {
        const int max_iter = 20;
        const int min_iter = 1;
        bool is_converged = false;
        // solve using Newton
        for (int iter = 0; iter < max_iter; ++iter) {
            assembleEquations();
            if (iter >= min_iter && converged(iter)) {
                is_converged = true;
                break;
            }
            solveAndUpdate();
        }
        if (!is_converged) {
            const auto msg =
                fmt::format(fmt::runtime("Temperature model (TEMP): Newton did not converge after {} iterations. \n"
                                         "The Simulator will continue to the next step with an unconverged solution."),
                            max_iter);
            OpmLog::debug(msg);
        }
    }

    void solveAndUpdate()
    {
        const unsigned int numCells = simulator_.model().numTotalDof();
        EnergyVector dx(numCells);
        bool conv = this->linearSolve_(*this->energyMatrix_, dx, this->energyVector_);
        if (!conv) {
            if (simulator_.gridView().comm().rank() == 0) {
                OpmLog::warning("Temp model: Linear solver did not converge. Temperature values not updated.");
            }
        }
        else {
            for (unsigned globI = 0; globI < numCells; ++globI) {
                this->temperature_[globI] -= std::clamp(dx[globI][0], -this->maxTempChange_, this->maxTempChange_);
                intQuants_[globI].updateTemperature_(simulator_.problem(), globI, /*timeIdx*/ 0);
                intQuants_[globI].updateEnergyQuantities_(simulator_.problem(), globI, /*timeIdx*/ 0);
            }
        }
    }

    bool converged(const int iter)
    {
        Scalar dt = simulator_.timeStepSize();
        Scalar maxNorm = 0.0;
        Scalar sumNorm = 0.0;
        const auto tolerance_cnv_energy_strict = Parameters::Get<Parameters::ToleranceCnvEnergy<Scalar>>();
        const auto& elemMapper = simulator_.model().elementMapper();
        const IsNumericalAquiferCell isNumericalAquiferCell(simulator_.gridView().grid());
        Scalar sum_pv = 0.0;
        Scalar sum_pv_not_converged = 0.0;
        for (const auto& elem : elements(simulator_.gridView(), Dune::Partitions::interior)) {
            unsigned globI = elemMapper.index(elem);
            const auto pvValue =  simulator_.problem().referencePorosity(globI, /*timeIdx=*/0)
            *  simulator_.model().dofTotalVolume(globI);

            const Scalar scaled_norm = dt * std::abs(this->energyVector_[globI])/ pvValue;
            maxNorm = max(maxNorm, scaled_norm);
            sumNorm += scaled_norm;
            if (!isNumericalAquiferCell(elem)) {
                if (scaled_norm > tolerance_cnv_energy_strict) {
                    sum_pv_not_converged += pvValue;
                }
                sum_pv += pvValue;
            }
        }
        maxNorm = simulator_.gridView().comm().max(maxNorm);
        sumNorm = simulator_.gridView().comm().sum(sumNorm);
        sum_pv = simulator_.gridView().comm().sum(sum_pv);
        sumNorm /= sum_pv;

        // Use relaxed tolerance if the fraction of unconverged cells porevolume is less than relaxed_max_pv_fraction
        sum_pv_not_converged = simulator_.gridView().comm().sum(sum_pv_not_converged);
        Scalar relaxed_max_pv_fraction = Parameters::Get<Parameters::RelaxedMaxPvFraction<Scalar>>();
        const bool relax = (sum_pv_not_converged / sum_pv) <  relaxed_max_pv_fraction;
        const auto tolerance_energy_balance = relax? Parameters::Get<Parameters::ToleranceEnergyBalanceRelaxed<Scalar>>():
                                            Parameters::Get<Parameters::ToleranceEnergyBalance<Scalar>>();
        const bool tolerance_cnv_energy = relax? Parameters::Get<Parameters::ToleranceCnvEnergyRelaxed<Scalar>>():
                                            tolerance_cnv_energy_strict;

        const auto msg = fmt::format(fmt::runtime("Temperature model (TEMP): Newton iter {}: "
                                     "CNV(E): {:.1e}, EB: {:.1e}"),
                                     iter, maxNorm, sumNorm);
        OpmLog::debug(msg);
        if (maxNorm < tolerance_cnv_energy && sumNorm < tolerance_energy_balance) {
            const auto msg2 =
                fmt::format(fmt::runtime("Temperature model (TEMP): Newton converged after {} iterations"),
                                         iter);
            OpmLog::debug(msg2);
            return true;
        }
        return false;
    }

    template< class LhsEval>
    void computeStorageTerm(unsigned globI, LhsEval& storage)
    {
        const auto& intQuants = intQuants_[globI];
        const auto& poro = decay<LhsEval>(intQuants.porosity());
        // accumulate the internal energy of the fluids
        const auto& fs = intQuants.fluidState();
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            const auto& u = decay<LhsEval>(fs.internalEnergy(phaseIdx));
            const auto& S = decay<LhsEval>(fs.saturation(phaseIdx));
            const auto& rho = decay<LhsEval>(fs.density(phaseIdx));

            storage += poro*S*u*rho;
        }

        // add the internal energy of the rock
        const Scalar rockFraction = intQuants.rockFraction();
        const auto& uRock = decay<LhsEval>(intQuants.rockInternalEnergy());
        storage += rockFraction*uRock;
        storage*= getPropValue<TypeTag, Properties::BlackOilEnergyScalingFactor>();
    }

    template < class ResidualNBInfo>
    void computeFluxTerm(unsigned globI, unsigned globJ,
                         const ResidualNBInfo& res_nbinfo,
                         Evaluation& flux)
    {
        const IntensiveQuantities& intQuantsIn = intQuants_[globI];
        const IntensiveQuantities& intQuantsEx = intQuants_[globJ];
        RateVector tmp(0.0); //not used
        RateVector darcyFlux(0.0);
        LocalResidual::computeFlux(tmp, darcyFlux, globI, globJ, intQuantsIn, intQuantsEx,
                                   res_nbinfo, simulator_.problem().moduleParams());
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            const unsigned activeCompIdx =
                FluidSystem::canonicalToActiveCompIdx(FluidSystem::solventComponentIndex(phaseIdx));

            bool inIsUp = darcyFlux[activeCompIdx] > 0;
            const IntensiveQuantities& up = inIsUp ? intQuantsIn : intQuantsEx;
            const auto& fs = up.fluidState();
            if (inIsUp) {
                flux += fs.enthalpy(phaseIdx)
                    * fs.density(phaseIdx)
                    * darcyFlux[activeCompIdx];
            }
            else {
                flux += getValue(fs.enthalpy(phaseIdx))
                    * getValue(fs.density(phaseIdx))
                    * getValue(darcyFlux[activeCompIdx]);
            }
        }
        flux *= getPropValue<TypeTag, Properties::BlackOilEnergyScalingFactor>();
    }

    template < class ResidualNBInfo>
    void computeHeatFluxTerm(unsigned globI, unsigned globJ,
                             const ResidualNBInfo& res_nbinfo,
                             Evaluation& heatFlux)
    {
        const IntensiveQuantities& intQuantsIn = intQuants_[globI];
        const IntensiveQuantities& intQuantsEx = intQuants_[globJ];
        const Scalar inAlpha = simulator_.problem().thermalHalfTransmissibility(globI, globJ);
        const Scalar outAlpha = simulator_.problem().thermalHalfTransmissibility(globJ, globI);
        short interiorDofIdx = 0; // NB
        short exteriorDofIdx = 1; // NB
        EnergyModule::ExtensiveQuantities::updateEnergy(heatFlux,
                                                        interiorDofIdx, // focusDofIndex,
                                                        interiorDofIdx,
                                                        exteriorDofIdx,
                                                        intQuantsIn,
                                                        intQuantsEx,
                                                        intQuantsIn.fluidState(),
                                                        intQuantsEx.fluidState(),
                                                        inAlpha,
                                                        outAlpha,
                                                        res_nbinfo.faceArea);
        heatFlux *= getPropValue<TypeTag, Properties::BlackOilEnergyScalingFactor>()*res_nbinfo.faceArea;
    }

    void assembleEquations()
    {
        this->energyVector_ = 0.0;
        (*this->energyMatrix_) = 0.0;
        Scalar dt = simulator_.timeStepSize();
        const unsigned int numCells = simulator_.model().numTotalDof();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (unsigned globI = 0; globI < numCells; ++globI) {
            Scalar volume = simulator_.model().dofTotalVolume(globI);
            Scalar storefac = volume / dt;
            Evaluation storage = 0.0;
            computeStorageTerm(globI, storage);
            this->energyVector_[globI] += storefac * ( getValue(storage) - storage1_[globI] );
            (*this->energyMatrix_)[globI][globI][0][0] += storefac * storage.derivative(Indices::temperatureIdx);
        }

        const auto& neighborInfo = simulator_.model().linearizer().getNeighborInfo();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (unsigned globI = 0; globI < numCells; ++globI) {
            const auto& nbInfos = neighborInfo[globI];
            for (const auto& nbInfo : nbInfos) {
                unsigned globJ = nbInfo.neighbor;
                assert(globJ != globI);

                // compute convective flux
                Evaluation flux = 0.0;
                computeFluxTerm(globI, globJ, nbInfo.res_nbinfo, flux);
                this->energyVector_[globI] += getValue(flux);
                (*this->energyMatrix_)[globI][globI][0][0] += flux.derivative(Indices::temperatureIdx);
                (*this->energyMatrix_)[globJ][globI][0][0] -= flux.derivative(Indices::temperatureIdx);

                // compute conductive flux
                Evaluation heatFlux = 0.0;
                computeHeatFluxTerm(globI, globJ, nbInfo.res_nbinfo, heatFlux);
                this->energyVector_[globI] += getValue(heatFlux);
                (*this->energyMatrix_)[globI][globI][0][0] += heatFlux.derivative(Indices::temperatureIdx);
                (*this->energyMatrix_)[globJ][globI][0][0] -= heatFlux.derivative(Indices::temperatureIdx);
            }
        }

        // Well terms
        const auto& wellPtrs = simulator_.problem().wellModel().localNonshutWells();
        for (const auto& wellPtr : wellPtrs) {
            this->assembleEquationWell(*wellPtr);
        }

        const auto& problem = simulator_.problem();

        bool enableDriftCompensation = Parameters::Get<Parameters::EnableDriftCompensationTemp>();
        if (enableDriftCompensation) {
            for (unsigned globalDofIdx = 0; globalDofIdx < numCells; ++globalDofIdx) {
                auto dofDriftRate = problem.drift()[globalDofIdx]/dt;
                const auto& fs = intQuants_[globalDofIdx].fluidState();
                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                   const unsigned activeCompIdx =
                        FluidSystem::canonicalToActiveCompIdx(FluidSystem::solventComponentIndex(phaseIdx));
                   auto drift_hrate = dofDriftRate[activeCompIdx]*getValue(fs.enthalpy(phaseIdx)) * getValue(fs.density(phaseIdx)) /  getValue(fs.invB(phaseIdx));
                   this->energyVector_[globalDofIdx] -= drift_hrate*getPropValue<TypeTag, Properties::BlackOilEnergyScalingFactor>();
                }
            }
        }

        if (simulator_.gridView().comm().size() > 1) {
            // Set dirichlet conditions for overlapping cells
            // loop over precalculated overlap rows and columns
            for (const auto row : overlapRows_) {
                // Zero out row.
                (*this->energyMatrix_)[row] = 0.0;

                //diagonal block set to diag(1.0).
                (*this->energyMatrix_)[row][row][0][0] = 1.0;
            }
        }
    }

    template<class Well>
    void assembleEquationWell(const Well& well)
    {
        const auto& eclWell = well.wellEcl();
        std::size_t well_index = simulator_.problem().wellModel().wellState().index(well.name()).value();
        const auto& ws = simulator_.problem().wellModel().wellState().well(well_index);
        this->energy_rates_[well_index] = 0.0;
        for (std::size_t i = 0; i < ws.perf_data.size(); ++i) {
            const auto globI = ws.perf_data.cell_index[i];
            auto fs = intQuants_[globI].fluidState(); //copy to make it possible to change the temp in the injector
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                Evaluation rate = well.volumetricSurfaceRateForConnection(globI, phaseIdx);
                if (rate > 0 && eclWell.isInjector()) {
                    fs.setTemperature(eclWell.inj_temperature());
                    const auto& rho = FluidSystem::density(fs, phaseIdx, fs.pvtRegionIndex());
                    fs.setDensity(phaseIdx, rho);
                    const auto& h = FluidSystem::enthalpy(fs, phaseIdx, fs.pvtRegionIndex());
                    fs.setEnthalpy(phaseIdx, h);
                    rate *= getValue(fs.enthalpy(phaseIdx)) * getValue(fs.density(phaseIdx)) / getValue(fs.invB(phaseIdx));
                } else {
                    const Evaluation d =  1.0 - fs.Rv() * fs.Rs();
                    if (phaseIdx == gasPhaseIdx && d > 0) {
                        const auto& oilrate = well.volumetricSurfaceRateForConnection(globI, oilPhaseIdx);
                        rate -= oilrate * getValue(fs.Rs());
                        rate /= d;
                    }
                    if (phaseIdx == oilPhaseIdx && d > 0) {
                        const auto& gasrate = well.volumetricSurfaceRateForConnection(globI, gasPhaseIdx);
                        rate -= gasrate * getValue(fs.Rv());
                        rate /= d;
                    }
                    rate *= fs.enthalpy(phaseIdx) * getValue(fs.density(phaseIdx)) / getValue(fs.invB(phaseIdx));
                }
                this->energy_rates_[well_index] += getValue(rate);
                rate *= getPropValue<TypeTag, Properties::BlackOilEnergyScalingFactor>();
                this->energyVector_[globI] -= getValue(rate);
                (*this->energyMatrix_)[globI][globI][0][0] -= rate.derivative(Indices::temperatureIdx);
            }
        }
    }

    const Simulator& simulator_;
    EnergyVector storage1_;
    std::vector<IntensiveQuantities> intQuants_;
    std::vector<int> overlapRows_;
    std::vector<int> interiorRows_;
};

// need for the old linearizer
template <class TypeTag>
class TemperatureModel<TypeTag, false>
{
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

public:
    TemperatureModel(Simulator&)
    { }

    /*!
     * \brief This method writes the complete state of all temperature
     *        to the hard disk.
     */
    template <class Restarter>
    void serialize(Restarter&)
    { /* not implemented */ }

    /*!
     * \brief This method restores the complete state of the temperature
     *        from disk.
     *
     * It is the inverse of the serialize() method.
     */
    template <class Restarter>
    void deserialize(Restarter&)
    { /* not implemented */ }

    void init() {}
    void beginTimeStep() {}
    const Scalar temperature(size_t /*globalIdx*/) const
    {
        return 273.15; // return 0C to make the compiler happy
    }
};

} // namespace Opm

#endif // OPM_TEMPERATURE_MODEL_HPP
