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
        this->doInit(rst, simulator_.model().numGridDof());

        if (!this->doTemp())
            return;

        storage0_.resize(simulator_.model().numGridDof());
    }

    void beginTimeStep()
    {
        if (!this->doTemp())
            return;

        updateStorageCache();
    }

    /*!
     * \brief Informs the temperature model that a time step has just been finished.
     */
    void endTimeStep()
    {
        if (!this->doTemp())
            return;

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
            //computeStorageTerm(globI, storage);
            storage0_[globI] = storage;
        }
        std::cout << "updateStorageCache" << std::endl;
    }

    void advanceTemperatureFields()
    {

        const unsigned int numCells = simulator_.model().numTotalDof();
        //#ifdef _OPENMP
        //#pragma omp parallel for
        //#endif
        Scalar dt = simulator_.timeStepSize();
        for (unsigned globI = 0; globI < numCells; ++globI) {
            Scalar volume = simulator_.model().dofTotalVolume(globI);
            Scalar storefac = volume / dt;
            Scalar storage = 0.0;
            //computeStorageTerm(globI, storage);
            //TemperatureEvaluation tmpStorage(getValue(storage));
            //tmpStorage.setDerivative(storage.derivative[Indices::energyIdx]);
            this->energyVector_[globI] = storefac * ( storage0_[globI] - getValue(storage) );
        }
        std::cout << "advanceTemperatureFields" << std::endl;
    }

    template< class LhsEval>
    void computeStorageTerm(unsigned globI, LhsEval& storage) {
        const auto& intQuants = simulator_.model().intensiveQuantities(globI, /*timeIdx*/ 0);
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
   

    const Simulator& simulator_;
    EnergyVector storage0_;

};

} // namespace Opm

#endif // OPM_TEMPERATURE_MODEL_HPP
