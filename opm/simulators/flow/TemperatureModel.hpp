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

#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/GenericTemperatureModel.hpp>
//#include <opm/simulators/utils/VectorVectorDataHandle.hpp>

#include <array>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace Opm::Properties {

template<class TypeTag, class MyTypeTag>
struct EnableTemperatureModel {
    using type = UndefinedProperty;
};

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup BlackOilSimulator
 *
 * \brief A class which handles sequential implicit solution of the energy equation as specified in by TEMP
 */
template <class TypeTag>
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

    using EnergyMatrix = typename BaseType::EnergyMatrix;
    using EnergyVector = typename BaseType::EnergyVector;

    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { numPhases = FluidSystem::numPhases };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };

public:
    TemperatureModel(Simulator& simulator)
        : BaseType(simulator.vanguard().gridView(),
                   simulator.vanguard().eclState(),
                   simulator.vanguard().cartesianIndexMapper(),
                   simulator.model().dofMapper())
        , simulator_(simulator)
    { }

    void init(bool rst)
    {
        const unsigned int numCells = simulator_.model().numTotalDof();

        this->doInit(rst, numCells);

        if (!this->doTemp())
            return;

        storage0_.resize(numCells);

        for (unsigned globI = 0; globI < numCells; ++globI) {
            this->temperature_[globI] = simulator_.problem().initialFluidState(globI).temperature(0);
        }
     
        intQuants_.resize(numCells);
    }

    void beginTimeStep()
    {
        if (!this->doTemp())
            return;

        // We copy the intensive quantities here to make it possible to update them
        const unsigned int numCells = simulator_.model().numTotalDof();
        for (unsigned globI = 0; globI < numCells; ++globI) {
            intQuants_[globI] = simulator_.model().intensiveQuantities(globI, /*timeIdx*/ 0);
        }
        updateStorageCache();
    }

    /*!
     * \brief Informs the temperature model that a time step has just been finished.
     */
    void endTimeStep()
    {
        if (!this->doTemp())
            return;

        // We copy the intensive quantities here to make it possible to update them
        const unsigned int numCells = simulator_.model().numTotalDof();
        for (unsigned globI = 0; globI < numCells; ++globI) {
            intQuants_[globI] = simulator_.model().intensiveQuantities(globI, /*timeIdx*/ 0);
        }
        advanceTemperatureFields();
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

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(static_cast<BaseType&>(*this));
    }

protected:

    void updateStorageCache()
    {
        const unsigned int numCells = simulator_.model().numTotalDof();
        //#ifdef _OPENMP
        //#pragma omp parallel for
        //#endif
        for (unsigned globI = 0; globI < numCells; ++globI) {
            Scalar storage = 0.0;
            computeStorageTerm(globI, storage);
            storage0_[globI] = storage;
        }
        std::cout << "updateStorageCache" << std::endl;
    }

    void advanceTemperatureFields()
    {


        for (int iter = 0; iter < 20; ++iter) {
        this->energyVector_ = 0.0;
        (*this->energyMatrix_) = 0.0;

        Scalar dt = simulator_.timeStepSize();
        const unsigned int numCells = simulator_.model().numTotalDof();
        //#ifdef _OPENMP
        //#pragma omp parallel for
        //#endif
        for (unsigned globI = 0; globI < numCells; ++globI) {
            Scalar volume = simulator_.model().dofTotalVolume(globI);
            Scalar storefac = volume / dt;
            Evaluation storage = 0.0;
            computeStorageTerm(globI, storage); 
            this->energyVector_[globI] += storefac * ( getValue(storage) - storage0_[globI] );
            if (globI == 0) {
                std::cout << "storage " << getValue(storage) << " " << storage0_[globI] << " " << storefac * ( getValue(storage) - storage0_[globI] ) << " " <<  storefac * storage.derivative(Indices::temperatureIdx) << " " << dt/(3600*24) << std::endl;
            }
            (*this->energyMatrix_)[globI][globI][0][0] += storefac * storage.derivative(Indices::temperatureIdx); 
        }

        const auto& neighborInfo = simulator_.model().linearizer().getNeighborInfo();
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
        for (unsigned globI = 0; globI < numCells; ++globI) {
            RateVector tmp(0.0);
            RateVector darcyFlux(0.0);
            const auto& nbInfos = neighborInfo[globI];
            const IntensiveQuantities& intQuantsIn = intQuants_[globI];
            for (const auto& nbInfo : nbInfos) {
                unsigned globJ = nbInfo.neighbor;
                assert(globJ != globI);
                tmp = 0.0;
                darcyFlux = 0.0;
                const IntensiveQuantities& intQuantsEx = intQuants_[globJ];
                LocalResidual::computeFlux(tmp, darcyFlux, globI, globJ, intQuantsIn, intQuantsEx, nbInfo.res_nbinfo, simulator_.problem().moduleParams());
                //adres *= nbInfo.res_nbinfo.faceArea;

                Evaluation flux = 0.0;
                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                    if (!FluidSystem::phaseIsActive(phaseIdx))
                        continue;
                    unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                    bool inIsUp = darcyFlux[activeCompIdx] > 0;
                    const IntensiveQuantities& up = inIsUp ? intQuantsIn : intQuantsEx;
                    const auto& fs = up.fluidState();
                    if (inIsUp) {
                        flux += fs.enthalpy(phaseIdx)
                         * fs.density(phaseIdx)
                         * darcyFlux[activeCompIdx];
                    } else {
                        flux += getValue(fs.enthalpy(phaseIdx))
                         * getValue(fs.density(phaseIdx))
                         * getValue(darcyFlux[activeCompIdx]);
                    }
                }
                flux*= getPropValue<TypeTag, Properties::BlackOilEnergyScalingFactor>();
            
                this->energyVector_[globI] += getValue(flux);
                if (globI == 0) {
                    std::cout << "flux: " << globJ << " " << flux << std::endl;
                }

                (*this->energyMatrix_)[globI][globI][0][0] += flux.derivative(Indices::temperatureIdx);
                (*this->energyMatrix_)[globJ][globI][0][0] -= flux.derivative(Indices::temperatureIdx);

            }
        }

        // Well terms
        const auto& wellPtrs = simulator_.problem().wellModel().localNonshutWells();
        for (const auto& wellPtr : wellPtrs) {
            this->assembleEquationWell(*wellPtr);
        }
        Scalar maxNorm = 0.0;
        Scalar sumNorm = 0.0;
        for (unsigned globI = 0; globI < numCells; ++globI) {
            maxNorm = max(maxNorm, std::abs(this->energyVector_[globI]));
            sumNorm += std::abs(this->energyVector_[globI]);
        }

        if (maxNorm < 1e-3 || sumNorm/numCells < 1e-7) {
            std::cout << "Newton converged" << std::endl;
            break;
        } 
        std::cout << "iter continue: " << iter << " " << maxNorm << " " << sumNorm/numCells << std::endl;

        EnergyVector dx(numCells);
        bool converged = this->linearSolve_(*this->energyMatrix_, dx, this->energyVector_);
        if (!converged) {
            OpmLog::warning("### Temp model: Linear solver did not converge. ###");
        }
        for (unsigned globI = 0; globI < numCells; ++globI) {
            if (globI == 0) {
                std::cout << globI << " " << this->temperature_[globI] << " " << dx[globI] << " " << this->energyVector_[globI] << std::endl;
            }
            this->temperature_[globI] -= std::clamp(dx[globI][0], Scalar(-10.0), Scalar(10.0));

            intQuants_[globI].updateTemperature_(simulator_.problem(), globI, /*timeIdx*/ 0);
            intQuants_[globI].updateEnergyQuantities_(simulator_.problem(), globI, /*timeIdx*/ 0);
        }
        }
        std::cout << "advanceTemperatureFields" << std::endl;
    }

    template< class LhsEval>
    void computeStorageTerm(unsigned globI, LhsEval& storage) {
        const auto& intQuants = intQuants_[globI];
        const auto& poro = decay<LhsEval>(intQuants.porosity());
        // accumulate the internal energy of the fluids
        const auto& fs = intQuants.fluidState();
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            const auto& u = decay<LhsEval>(fs.internalEnergy(phaseIdx));
            const auto& S = decay<LhsEval>(fs.saturation(phaseIdx));
            const auto& rho = decay<LhsEval>(fs.density(phaseIdx));

            storage += poro*S*u*rho;
        }

        // add the internal energy of the rock
        Scalar rockFraction = intQuants.rockFraction();
        const auto& uRock = decay<LhsEval>(intQuants.rockInternalEnergy());
        storage += rockFraction*uRock;
        storage*= getPropValue<TypeTag, Properties::BlackOilEnergyScalingFactor>();
    }

    template<class Well>
    void assembleEquationWell(const Well& well)
    {
        const auto& eclWell = well.wellEcl();
        Scalar perf_temp = 200;
        Scalar dt = simulator_.timeStepSize();
        std::size_t well_index = simulator_.problem().wellModel().wellState().index(well.name()).value();
        const auto& ws = simulator_.problem().wellModel().wellState().well(well_index);
        for (std::size_t i = 0; i < ws.perf_data.size(); ++i) {
            const auto globI = ws.perf_data.cell_index[i];
            auto fs = intQuants_[globI].fluidState(); //copy to make it possible to change the temp in the injector
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                //Scalar rate = ws.perf_data.rates[i]; //well.volumetricSurfaceRateForConnection(globI, phaseIdx);
                Scalar rate = well.volumetricSurfaceRateForConnection(globI, phaseIdx);

                if (phaseIdx == gasPhaseIdx) { // assumes rv == 0
                    rate -= well.volumetricSurfaceRateForConnection(globI, oilPhaseIdx)*getValue(fs.Rs());
                }

                if (rate > 0) {
                    fs.setTemperature(eclWell.inj_temperature());
                    const auto& rho = FluidSystem::density(fs, phaseIdx, fs.pvtRegionIndex());
                    fs.setDensity(phaseIdx, rho);
                    const auto& h = FluidSystem::enthalpy(fs, phaseIdx, fs.pvtRegionIndex());
                    fs.setEnthalpy(phaseIdx, h);
                    rate *= getValue(fs.enthalpy(phaseIdx)) * getValue(fs.density(phaseIdx)) / getValue(fs.invB(phaseIdx));
                } else {
                    rate *= getValue(fs.enthalpy(phaseIdx)) * getValue(fs.density(phaseIdx)) / getValue(fs.invB(phaseIdx));
                }
                std::cout << "temp well " << rate << std::endl;
                this->energyVector_[globI] -= rate * getPropValue<TypeTag, Properties::BlackOilEnergyScalingFactor>();
            }
        }
    }
   

    Simulator& simulator_;
    EnergyVector storage0_;
    std::vector<IntensiveQuantities> intQuants_;

};

} // namespace Opm

#endif // OPM_TEMPERATURE_MODEL_HPP
