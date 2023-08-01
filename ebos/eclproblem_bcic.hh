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
 * \copydoc Opm::EclProblem
 */
#ifndef OPM_ECL_PROBLEM_BCIC_HH
#define OPM_ECL_PROBLEM_BCIC_HH

#include <ebos/eclequilinitializer.hh>
#include <ebos/ecloutputblackoilmodule.hh>
#include <ebos/eclsolutioncontainers.hh>

#include <opm/common/TimingMacros.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/input/eclipse/Schedule/BCProp.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Opm {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Handling of boundary- and initial conditions for EclProblem.
 */
template <class TypeTag>
class EclProblemBCIC
{
public:
    using EclMaterialLawManager = typename GetProp<TypeTag, Properties::MaterialLaw>::EclMaterialLawManager;
    using MaterialLawParams = typename EclMaterialLawManager::MaterialLawParams;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using InitialFluidState = typename EclEquilInitializer<TypeTag>::ScalarFluidState;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Vanguard = GetPropType<TypeTag, Properties::Vanguard>;

    static constexpr bool enableBrine = getPropValue<TypeTag, Properties::EnableBrine>();
    static constexpr bool enableMICP = getPropValue<TypeTag, Properties::EnableMICP>();
    static constexpr bool enablePolymer = getPropValue<TypeTag, Properties::EnablePolymer>();
    static constexpr bool enablePolymerMW = getPropValue<TypeTag, Properties::EnablePolymerMW>();
    static constexpr bool enableSaltPrecipitation = getPropValue<TypeTag, Properties::EnableSaltPrecipitation>();
    static constexpr bool enableSolvent = getPropValue<TypeTag, Properties::EnableSolvent>();

    //! \brief Returns whether or not boundary conditions are trivial.
    bool nonTrivialBoundaryConditions() const
    {
        return nonTrivialBoundaryConditions_;
    }

    //! \brief Returns a const reference to initial fluid state for an element.
    const InitialFluidState& initialFluidState(const unsigned idx) const
    {
        return initialFluidStates_[idx];
    }

    //! \brief Reads boundary conditions from eclipse state.
    void readBoundaryConditions(const Vanguard& vanguard)
    {
        const auto& bcconfig = vanguard.eclState().getSimulationConfig().bcconfig();
        if (bcconfig.size() > 0) {
            nonTrivialBoundaryConditions_ = true;

            std::size_t numCartDof = vanguard.cartesianSize();
            unsigned numElems = vanguard.gridView().size(/*codim=*/0);
            std::vector<int> cartesianToCompressedElemIdx(numCartDof, -1);

            for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
                cartesianToCompressedElemIdx[vanguard.cartesianIndex(elemIdx)] = elemIdx;
            }

            bcindex_.resize(numElems, 0);
            auto loopAndApply = [&cartesianToCompressedElemIdx,
                                 &vanguard](const auto& bcface,
                                            auto apply)
            {
                for (int i = bcface.i1; i <= bcface.i2; ++i) {
                    for (int j = bcface.j1; j <= bcface.j2; ++j) {
                        for (int k = bcface.k1; k <= bcface.k2; ++k) {
                            std::array<int, 3> tmp = {i,j,k};
                            auto elemIdx = cartesianToCompressedElemIdx[vanguard.cartesianIndex(tmp)];
                            if (elemIdx >= 0) {
                                apply(elemIdx);
                            }
                        }
                    }
                }
            };
            for (const auto& bcface : bcconfig) {
                std::vector<int>& data = bcindex_(bcface.dir);
                const int index = bcface.index;
                loopAndApply(bcface,
                             [&data,index](int elemIdx)
                             { data[elemIdx] = index; });
            }
        }
    }

    //! \brief Reads the configured boundary conditions from eclipse state.
    void readInitialCondition(EclMaterialLawManager& materialLawManager,
                              [[maybe_unused]] std::vector<Scalar>& solventSaturation,
                              [[maybe_unused]] std::vector<Scalar>& solventRsw,
                              [[maybe_unused]] MICPSolutionContainer<Scalar>& micp,
                              [[maybe_unused]] PolymerSolutionContainer<Scalar>& polymer,
                              const Simulator& simulator,
                              const std::size_t numGridDof,
                              const std::function<int(int)> pvtRegionIndex)
    {
        const auto& eclState = simulator.vanguard().eclState();

        if (eclState.getInitConfig().hasEquil()) {
            this->readEquilInitialCondition_(materialLawManager,
                                             simulator,
                                             numGridDof);
        } else {
            this->readExplicitInitialCondition_(eclState.fieldProps(),
                                                materialLawManager,
                                                numGridDof,
                                                pvtRegionIndex);
        }

        if constexpr (enableSolvent) {
            this->readSolventInitialCondition_(solventSaturation,
                                               solventRsw,
                                               eclState.fieldProps(),
                                               numGridDof);
        }

        if constexpr (enableMICP) {
            micp.readInitialCondition(eclState.fieldProps(), numGridDof);
        }

        if constexpr (enablePolymer) {
            polymer.readInitialCondition(eclState.fieldProps(),
                                         enablePolymerMW, numGridDof);
        }
    }

    //! \brief Read restart solution from eclipse input.
    void readEclRestartSolution([[maybe_unused]] std::vector<Scalar>& solventSaturation,
                                [[maybe_unused]] std::vector<Scalar>& solventRsw,
                                [[maybe_unused]] MICPSolutionContainer<Scalar>& micp,
                                [[maybe_unused]] PolymerSolutionContainer<Scalar>& polymer,
                                const EclOutputBlackOilModule<TypeTag>& input,
                                const std::size_t numElems,
                                const std::function<int(int)> pvtRegionIndex)
    {
        initialFluidStates_.resize(numElems);
        for (std::size_t elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            initialFluidStates_[elemIdx].setPvtRegionIndex(pvtRegionIndex(elemIdx));
            input.assignToFluidState(initialFluidStates_[elemIdx], elemIdx);
            // Note: Function processRestartSaturations_() mutates the
            // 'ssol' argument--the value from the restart file--if solvent
            // is enabled.  Then, store the updated solvent saturation into
            // 'solventSaturation_'.  Otherwise, just pass a dummy value to
            // the function and discard the unchanged result.  Do not index
            // into 'solventSaturation_' unless solvent is enabled.
            auto ssol = enableSolvent
                ? input.getSolventSaturation(elemIdx)
                : Scalar(0);

            processRestartSaturations_(elemIdx, ssol);

            if constexpr (enableSolvent) {
                solventSaturation[elemIdx] = ssol;
                solventRsw[elemIdx] = input.getSolventRsw(elemIdx);
            }

            if constexpr (enablePolymer) {
                polymer.concentration[elemIdx] = input.getPolymerConcentration(elemIdx);
            }

            // if we need to restart for polymer molecular weight simulation, we need to add related here

            if constexpr (enableMICP) {
                 micp.microbialConcentration[elemIdx] = input.getMicrobialConcentration(elemIdx);
                 micp.oxygenConcentration[elemIdx] = input.getOxygenConcentration(elemIdx);
                 micp.ureaConcentration[elemIdx] = input.getUreaConcentration(elemIdx);
                 micp.biofilmConcentration[elemIdx] = input.getBiofilmConcentration(elemIdx);
                 micp.calciteConcentration[elemIdx] = input.getCalciteConcentration(elemIdx);
            }
        }
    }

    InitialFluidState boundaryFluidState(unsigned globalDofIdx,
                                         const int directionId,
                                         const int pvtRegionIdx,
                                         const BCProp& bcprop,
                                         const MaterialLawParams& matParams) const
    {
        OPM_TIMEBLOCK_LOCAL(boundaryFluidState);
        if (bcprop.size() > 0) {
            FaceDir::DirEnum dir = FaceDir::FromIntersectionIndex(directionId);

            // index == 0: no boundary conditions for this
            // global cell and direction
            if (bcindex_(dir)[globalDofIdx] == 0) {
                return initialFluidState(globalDofIdx);
            }

            const auto& bc = bcprop[bcindex_(dir)[globalDofIdx]];
            if (bc.bctype == BCType::DIRICHLET )
            {
                InitialFluidState fluidState;
                fluidState.setPvtRegionIndex(pvtRegionIdx);

                switch (bc.component) {
                    case BCComponent::OIL:
                        if (!FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx))
                            throw std::logic_error("oil is not active and you're trying to add oil BC");

                        fluidState.setSaturation(FluidSystem::oilPhaseIdx, 1.0);
                        break;
                    case BCComponent::GAS:
                        if (!FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx))
                            throw std::logic_error("gas is not active and you're trying to add gas BC");

                        fluidState.setSaturation(FluidSystem::gasPhaseIdx, 1.0);
                        break;
                        case BCComponent::WATER:
                        if (!FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx))
                            throw std::logic_error("water is not active and you're trying to add water BC");

                        fluidState.setSaturation(FluidSystem::waterPhaseIdx, 1.0);
                        break;
                    case BCComponent::SOLVENT:
                    case BCComponent::POLYMER:
                    case BCComponent::NONE:
                        throw std::logic_error("you need to specify a valid component (OIL, WATER or GAS) when DIRICHLET type is set in BC");
                        break;
                }
                int phaseIndex;
                if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                    phaseIndex = FluidSystem::oilPhaseIdx;
                }
                else if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                    phaseIndex = FluidSystem::gasPhaseIdx;
                }
                else {
                    phaseIndex = FluidSystem::waterPhaseIdx;
                }
                double pressure = initialFluidState(globalDofIdx).pressure(phaseIndex);
                const auto pressure_input = bc.pressure;
                if (pressure_input) {
                    pressure = *pressure_input;
                }

                std::array<Scalar, FluidSystem::numPhases> pc = {0};
                MaterialLaw::capillaryPressures(pc, matParams, fluidState);
                Valgrind::CheckDefined(pressure);
                Valgrind::CheckDefined(pc);
                for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                    if (!FluidSystem::phaseIsActive(phaseIdx))
                        continue;

                    if (Indices::oilEnabled)
                        fluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[FluidSystem::oilPhaseIdx]));
                    else if (Indices::gasEnabled)
                        fluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[FluidSystem::gasPhaseIdx]));
                    else if (Indices::waterEnabled)
                        //single (water) phase
                        fluidState.setPressure(phaseIdx, pressure);
                }

                double temperature = initialFluidState(globalDofIdx).temperature(phaseIndex);
                const auto temperature_input = bc.temperature;
                if (temperature_input) {
                    temperature = *temperature_input;
                }
                fluidState.setTemperature(temperature);

                if (FluidSystem::enableDissolvedGas()) {
                    fluidState.setRs(0.0);
                    fluidState.setRv(0.0);
                }
                if (FluidSystem::enableDissolvedGasInWater()) {
                    fluidState.setRsw(0.0);
                }
                if (FluidSystem::enableVaporizedWater())
                    fluidState.setRvw(0.0);

                for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                    if (!FluidSystem::phaseIsActive(phaseIdx))
                        continue;

                    const auto& b = FluidSystem::inverseFormationVolumeFactor(fluidState, phaseIdx, pvtRegionIdx);
                    fluidState.setInvB(phaseIdx, b);

                    const auto& rho = FluidSystem::density(fluidState, phaseIdx, pvtRegionIdx);
                    fluidState.setDensity(phaseIdx, rho);

                }
                fluidState.checkDefined();
                return fluidState;
            }
        }
        return initialFluidState(globalDofIdx);
    }

    //! \brief Calculate equilibrium boundary conditions.
    void readEquilInitialCondition_(EclMaterialLawManager& materialLawManager,
                                    const Simulator& simulator,
                                    const std::size_t numElems)
    {
        // initial condition corresponds to hydrostatic conditions.
        using EquilInitializer = EclEquilInitializer<TypeTag>;
        EquilInitializer equilInitializer(simulator, materialLawManager);

        initialFluidStates_.resize(numElems);
        for (std::size_t elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            auto& elemFluidState = initialFluidStates_[elemIdx];
            elemFluidState.assign(equilInitializer.initialFluidState(elemIdx));
        }
    }

    //! \brief Read explicitly specified initial conditions from eclipse state.
    void readExplicitInitialCondition_(const FieldPropsManager& fp,
                                       const EclMaterialLawManager& materialLawManager,
                                       const std::size_t numDof,
                                       std::function<int(int)> pvtRegionIndex)
    {
        bool has_swat     = fp.has_double("SWAT");
        bool has_sgas     = fp.has_double("SGAS");
        bool has_rs       = fp.has_double("RS");
        bool has_rv       = fp.has_double("RV");
        bool has_rvw       = fp.has_double("RVW");
        bool has_pressure = fp.has_double("PRESSURE");
        bool has_salt = fp.has_double("SALT");
        bool has_saltp = fp.has_double("SALTP");

        // make sure all required quantities are enables
        if (Indices::numPhases > 1) {
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && !has_swat) {
                throw std::runtime_error("The ECL input file requires the presence of the SWAT keyword if "
                                     "the water phase is active");
            }
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) &&
                FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && !has_sgas) {
                throw std::runtime_error("The ECL input file requires the presence of the SGAS keyword if "
                                     "the gas phase is active");
            }
        }
        if (!has_pressure) {
            throw std::runtime_error("The ECL input file requires the presence of the PRESSURE "
                                      "keyword if the model is initialized explicitly");
        }
        if (FluidSystem::enableDissolvedGas() && !has_rs) {
            throw std::runtime_error("The ECL input file requires the RS keyword to be present if"
                                     " dissolved gas is enabled");
        }
        if (FluidSystem::enableVaporizedOil() && !has_rv) {
            throw std::runtime_error("The ECL input file requires the RV keyword to be present if"
                                     " vaporized oil is enabled");
        }
        if (FluidSystem::enableVaporizedWater() && !has_rvw) {
            throw std::runtime_error("The ECL input file requires the RVW keyword to be present if"
                                     " vaporized water is enabled");
        }
        if (enableBrine && !has_salt) {
            throw std::runtime_error("The ECL input file requires the SALT keyword to be present if"
                                     " brine is enabled and the model is initialized explicitly");
        }
        if (enableSaltPrecipitation && !has_saltp) {
            throw std::runtime_error("The ECL input file requires the SALTP keyword to be present if"
                                     " salt precipitation is enabled and the model is initialized explicitly");
        }

        initialFluidStates_.resize(numDof);

        std::vector<double> waterSaturationData;
        std::vector<double> gasSaturationData;
        std::vector<double> pressureData;
        std::vector<double> rsData;
        std::vector<double> rvData;
        std::vector<double> rvwData;
        std::vector<double> tempiData;
        std::vector<double> saltData;
        std::vector<double> saltpData;

        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && Indices::numPhases > 1) {
            waterSaturationData = fp.get_double("SWAT");
        } else {
            waterSaturationData.resize(numDof);
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) &&
            FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            gasSaturationData = fp.get_double("SGAS");
        } else {
            gasSaturationData.resize(numDof);
        }

        pressureData = fp.get_double("PRESSURE");
        if (FluidSystem::enableDissolvedGas()) {
            rsData = fp.get_double("RS");
        }

        if (FluidSystem::enableVaporizedOil()) {
            rvData = fp.get_double("RV");
        }

        if (FluidSystem::enableVaporizedWater()) {
            rvwData = fp.get_double("RVW");
        }

        // initial reservoir temperature
        tempiData = fp.get_double("TEMPI");

        // initial salt concentration data
        if constexpr (enableBrine) {
            saltData = fp.get_double("SALT");
        }

         // initial precipitated salt saturation data
        if constexpr (enableSaltPrecipitation) {
            saltpData = fp.get_double("SALTP");
        }

        // calculate the initial fluid states
        for (std::size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            auto& dofFluidState = initialFluidStates_[dofIdx];

            dofFluidState.setPvtRegionIndex(pvtRegionIndex(dofIdx));

            //////
            // set temperature
            //////
            Scalar temperatureLoc = tempiData[dofIdx];
            if (!std::isfinite(temperatureLoc) || temperatureLoc <= 0) {
                temperatureLoc = FluidSystem::surfaceTemperature;
            }
            dofFluidState.setTemperature(temperatureLoc);

            //////
            // set salt concentration
            //////
            if constexpr (enableBrine) {
                dofFluidState.setSaltConcentration(saltData[dofIdx]);
            }

            //////
            // set precipitated salt saturation
            //////
            if constexpr (enableSaltPrecipitation) {
                dofFluidState.setSaltSaturation(saltpData[dofIdx]);
            }

            //////
            // set saturations
            //////
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                dofFluidState.setSaturation(FluidSystem::waterPhaseIdx,
                                            waterSaturationData[dofIdx]);
            }

            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                if (!FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                    dofFluidState.setSaturation(FluidSystem::gasPhaseIdx,
                                            1.0
                                            - waterSaturationData[dofIdx]);
                } else {
                    dofFluidState.setSaturation(FluidSystem::gasPhaseIdx,
                                                gasSaturationData[dofIdx]);
                }
            }
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                dofFluidState.setSaturation(FluidSystem::oilPhaseIdx,
                                            1.0
                                            - waterSaturationData[dofIdx]
                                            - gasSaturationData[dofIdx]);
            }

            //////
            // set phase pressures
            //////
            Scalar pressure = pressureData[dofIdx]; // oil pressure (or gas pressure for water-gas system or water pressure for single phase)

            // this assumes that capillary pressures only depend on the phase saturations
            // and possibly on temperature. (this is always the case for ECL problems.)
            std::array<Scalar, FluidSystem::numPhases> pc = {0};
            const auto& matParams = materialLawManager.materialLawParams(dofIdx);
            MaterialLaw::capillaryPressures(pc, matParams, dofFluidState);
            Valgrind::CheckDefined(pressure);
            Valgrind::CheckDefined(pc);
            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                if (Indices::oilEnabled) {
                    dofFluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[FluidSystem::oilPhaseIdx]));
                } else if (Indices::gasEnabled) {
                    dofFluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[FluidSystem::gasPhaseIdx]));
                } else if (Indices::waterEnabled) {
                    //single (water) phase
                    dofFluidState.setPressure(phaseIdx, pressure);
                }
            }

            if (FluidSystem::enableDissolvedGas()) {
                dofFluidState.setRs(rsData[dofIdx]);
            } else if (Indices::gasEnabled && Indices::oilEnabled) {
                dofFluidState.setRs(0.0);
            }

            if (FluidSystem::enableVaporizedOil()) {
                dofFluidState.setRv(rvData[dofIdx]);
            } else if (Indices::gasEnabled && Indices::oilEnabled) {
                dofFluidState.setRv(0.0);
            }

            if (FluidSystem::enableVaporizedWater()) {
                dofFluidState.setRvw(rvwData[dofIdx]);
            }

            //////
            // set invB_
            //////
            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                const auto& b = FluidSystem::inverseFormationVolumeFactor(dofFluidState,
                                                                          phaseIdx,
                                                                          pvtRegionIndex(dofIdx));
                dofFluidState.setInvB(phaseIdx, b);

                const auto& rho = FluidSystem::density(dofFluidState, phaseIdx,
                                                       pvtRegionIndex(dofIdx));
                dofFluidState.setDensity(phaseIdx, rho);
            }
        }
    }

    //! \brief Read solvent initial condition from eclipse state.
    void readSolventInitialCondition_(std::vector<Scalar>& solventSaturation,
                                      std::vector<Scalar>& solventRsw,
                                      const FieldPropsManager& fp,
                                      const unsigned numDof)
    {
        if (fp.has_double("SSOL")) {
            solventSaturation = fp.get_double("SSOL");
        } else {
            solventSaturation.resize(numDof, 0.0);
        }
        solventRsw.resize(numDof, 0.0);
    }

    //! \brief Clamp small saturations due to single precision input.
    void processRestartSaturations_(const std::size_t elemIdx,
                                    Scalar& solventSaturation)
    {
        auto& elemFluidState = initialFluidStates_[elemIdx];
        // each phase needs to be above certain value to be claimed to be existing
        // this is used to recover some RESTART running with the defaulted single-precision format
        constexpr Scalar smallSaturationTolerance = 1.e-6;
        Scalar sumSaturation = 0.0;
        for (std::size_t phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
            if (FluidSystem::phaseIsActive(phaseIdx)) {
                if (elemFluidState.saturation(phaseIdx) < smallSaturationTolerance) {
                    elemFluidState.setSaturation(phaseIdx, 0.0);
                }

                sumSaturation += elemFluidState.saturation(phaseIdx);
            }
        }

        if constexpr (enableSolvent) {
            if (solventSaturation < smallSaturationTolerance) {
                solventSaturation = 0.0;
            }

           sumSaturation += solventSaturation;
        }

        assert(sumSaturation > 0.0);

        for (std::size_t phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
            if (FluidSystem::phaseIsActive(phaseIdx)) {
                const Scalar saturation = elemFluidState.saturation(phaseIdx) / sumSaturation;
                elemFluidState.setSaturation(phaseIdx, saturation);
            }
        }

        if constexpr (enableSolvent) {
            solventSaturation = solventSaturation / sumSaturation;
        }
    }

    //! \brief Container for boundary conditions.
    template<class T>
    struct BCData
    {
        std::array<std::vector<T>,6> data; //!< One vector per FaceDir::DirEnum entry

        //! \brief Resize vectors.
        void resize(std::size_t size, T defVal)
        {
            for (auto& d : data) {
                d.resize(size, defVal);
            }
        }

        //! \brief Returns a const reference to data for a direction.
        const std::vector<T>& operator()(FaceDir::DirEnum dir) const
        {
            if (dir == FaceDir::DirEnum::Unknown) {
                throw std::runtime_error("Tried to access BC data for the 'Unknown' direction");
            }
            int idx = 0;
            int div = static_cast<int>(dir);
            while ((div /= 2) >= 1) {
              ++idx;
            }
            assert(idx >= 0 && idx <= 5);
            return data[idx];
        }

        //! \brief Returns a reference to data for a direction.
        std::vector<T>& operator()(FaceDir::DirEnum dir)
        {
            return const_cast<std::vector<T>&>(std::as_const(*this)(dir));
        }
    };

    BCData<int> bcindex_; //!< Indices for boundary conditions
    bool nonTrivialBoundaryConditions_ = false; //!< Whether or not non-trivial boundary conditions are used
    std::vector<InitialFluidState> initialFluidStates_; //!< Vector of initial fluid states for elements
};

} // namespace Opm

#endif
