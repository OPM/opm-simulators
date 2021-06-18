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
 * \copydoc Opm::MultiPhaseBaseModel
 */
#ifndef EWOMS_MULTI_PHASE_BASE_MODEL_HH
#define EWOMS_MULTI_PHASE_BASE_MODEL_HH

#include <opm/material/densead/Math.hpp>

#include "multiphasebaseproperties.hh"
#include "multiphasebaseproblem.hh"
#include "multiphasebaseextensivequantities.hh"

#include <opm/models/common/flux.hh>
#include <opm/models/discretization/vcfv/vcfvdiscretization.hh>
#include <opm/models/io/vtkmultiphasemodule.hh>
#include <opm/models/io/vtktemperaturemodule.hh>

#include <opm/material/fluidmatrixinteractions/NullMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/thermal/NullThermalConductionLaw.hpp>
#include <opm/material/thermal/NullSolidEnergyLaw.hpp>
#include <opm/material/common/Unused.hpp>

namespace Opm {
template <class TypeTag>
class MultiPhaseBaseModel;
}

namespace Opm::Properties {

//! The generic type tag for problems using the immiscible multi-phase model
// Create new type tags
namespace TTag {
struct MultiPhaseBaseModel { using InheritsFrom = std::tuple<VtkTemperature, VtkMultiPhase>; };
} // end namespace TTag

//! Specify the splices of the MultiPhaseBaseModel type tag
template<class TypeTag>
struct Splices<TypeTag, TTag::MultiPhaseBaseModel>
{
    using type = std::tuple<GetSplicePropType<TypeTag, TTag::MultiPhaseBaseModel, Properties::SpatialDiscretizationSplice>>;
};

//! Set the default spatial discretization
//!
//! We use a vertex centered finite volume method by default
template<class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::MultiPhaseBaseModel> { using type = TTag::VcfvDiscretization; };

//! set the number of equations to the number of phases
template<class TypeTag>
struct NumEq<TypeTag, TTag::MultiPhaseBaseModel> { static constexpr int value = GetPropType<TypeTag, Properties::Indices>::numEq; };
//! The number of phases is determined by the fluid system
template<class TypeTag>
struct NumPhases<TypeTag, TTag::MultiPhaseBaseModel> { static constexpr int value = GetPropType<TypeTag, Properties::FluidSystem>::numPhases; };
//! Number of chemical species in the system
template<class TypeTag>
struct NumComponents<TypeTag, TTag::MultiPhaseBaseModel> { static constexpr int value = GetPropType<TypeTag, Properties::FluidSystem>::numComponents; };

//! The type of the base base class for actual problems
template<class TypeTag>
struct BaseProblem<TypeTag, TTag::MultiPhaseBaseModel> { using type = MultiPhaseBaseProblem<TypeTag>; };

//! By default, use the Darcy relation to determine the phase velocity
template<class TypeTag>
struct FluxModule<TypeTag, TTag::MultiPhaseBaseModel> { using type = DarcyFluxModule<TypeTag>; };

/*!
 * \brief Set the material law to the null law by default.
 */
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::MultiPhaseBaseModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Traits = NullMaterialTraits<Scalar, FluidSystem::numPhases>;

public:
    using type = NullMaterial<Traits>;
};

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
template<class TypeTag>
struct MaterialLawParams<TypeTag, TTag::MultiPhaseBaseModel>
{ using type = typename GetPropType<TypeTag, Properties::MaterialLaw>::Params; };

//! set the energy storage law for the solid to the one which assumes zero heat capacity
//! by default
template<class TypeTag>
struct SolidEnergyLaw<TypeTag, TTag::MultiPhaseBaseModel>
{ using type = NullSolidEnergyLaw<GetPropType<TypeTag, Properties::Scalar>>; };

//! extract the type of the parameter objects for the solid energy storage law from the
//! law itself
template<class TypeTag>
struct SolidEnergyLawParams<TypeTag, TTag::MultiPhaseBaseModel>
{ using type = typename GetPropType<TypeTag, Properties::SolidEnergyLaw>::Params; };

//! set the thermal conduction law to a dummy one by default
template<class TypeTag>
struct ThermalConductionLaw<TypeTag, TTag::MultiPhaseBaseModel>
{ using type = NullThermalConductionLaw<GetPropType<TypeTag, Properties::Scalar>>; };

//! extract the type of the parameter objects for the thermal conduction law from the law
//! itself
template<class TypeTag>
struct ThermalConductionLawParams<TypeTag, TTag::MultiPhaseBaseModel>
{ using type = typename GetPropType<TypeTag, Properties::ThermalConductionLaw>::Params; };

//! disable gravity by default
template<class TypeTag>
struct EnableGravity<TypeTag, TTag::MultiPhaseBaseModel> { static constexpr bool value = false; };


} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup MultiPhaseBaseModel
 * \brief A base class for fully-implicit multi-phase porous-media flow models
 *        which assume multiple fluid phases.
 */
template <class TypeTag>
class MultiPhaseBaseModel : public GetPropType<TypeTag, Properties::Discretization>
{
    using ParentType = GetPropType<TypeTag, Properties::Discretization>;
    using Implementation = GetPropType<TypeTag, Properties::Model>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using ThreadManager = GetPropType<TypeTag, Properties::ThreadManager>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

    using ElementIterator = typename GridView::template Codim<0>::Iterator;
    using Element = typename GridView::template Codim<0>::Entity;

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = FluidSystem::numComponents };

public:
    MultiPhaseBaseModel(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Register all run-time parameters for the immiscible model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        // register runtime parameters of the VTK output modules
        VtkMultiPhaseModule<TypeTag>::registerParameters();
        VtkTemperatureModule<TypeTag>::registerParameters();
    }

    /*!
     * \brief Returns true iff a fluid phase is used by the model.
     *
     * \param phaseIdx The index of the fluid phase in question
     */
    bool phaseIsConsidered(unsigned phaseIdx OPM_UNUSED) const
    { return true; }

    /*!
     * \brief Compute the total storage inside one phase of all
     *        conservation quantities.
     *
     * \copydetails Doxygen::storageParam
     * \copydetails Doxygen::phaseIdxParam
     */
    void globalPhaseStorage(EqVector& storage, unsigned phaseIdx)
    {
        assert(phaseIdx < numPhases);

        storage = 0;

        ThreadedEntityIterator<GridView, /*codim=*/0> threadedElemIt(this->gridView());
        std::mutex mutex;
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            // Attention: the variables below are thread specific and thus cannot be
            // moved in front of the #pragma!
            unsigned threadId = ThreadManager::threadId();
            ElementContext elemCtx(this->simulator_);
            ElementIterator elemIt = threadedElemIt.beginParallel();
            EqVector tmp;

            for (; !threadedElemIt.isFinished(elemIt); elemIt = threadedElemIt.increment()) {
                const Element& elem = *elemIt;
                if (elem.partitionType() != Dune::InteriorEntity)
                    continue; // ignore ghost and overlap elements

                elemCtx.updateStencil(elem);
                elemCtx.updateIntensiveQuantities(/*timeIdx=*/0);

                const auto& stencil = elemCtx.stencil(/*timeIdx=*/0);

                for (unsigned dofIdx = 0; dofIdx < elemCtx.numDof(/*timeIdx=*/0); ++dofIdx) {
                    const auto& scv = stencil.subControlVolume(dofIdx);
                    const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);

                    tmp = 0;
                    this->localResidual(threadId).addPhaseStorage(tmp,
                                                                  elemCtx,
                                                                  dofIdx,
                                                                  /*timeIdx=*/0,
                                                                  phaseIdx);
                    tmp *= scv.volume()*intQuants.extrusionFactor();

                    mutex.lock();
                    storage += tmp;
                    mutex.unlock();
                }
            }
        }

        storage = this->gridView_.comm().sum(storage);
    }

    void registerOutputModules_()
    {
        ParentType::registerOutputModules_();

        // add the VTK output modules which make sense for all multi-phase models
        this->addOutputModule(new VtkMultiPhaseModule<TypeTag>(this->simulator_));
        this->addOutputModule(new VtkTemperatureModule<TypeTag>(this->simulator_));
    }

private:
    const Implementation& asImp_() const
    { return *static_cast<const Implementation *>(this); }
};
} // namespace Opm

#endif
