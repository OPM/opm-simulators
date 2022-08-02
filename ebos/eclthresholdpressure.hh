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
 * \copydoc Opm::EclThresholdPressure
 */
#ifndef EWOMS_ECL_THRESHOLD_PRESSURE_HH
#define EWOMS_ECL_THRESHOLD_PRESSURE_HH

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/common/multiphasebaseproperties.hh>
#include <ebos/eclgenericthresholdpressure.hh>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>

#include <algorithm>
#include <vector>

#include <opm/models/discretization/common/smallelementcontext.hh>
namespace Opm {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief This class calculates the threshold pressure for grid faces according to the
 *        Eclipse Reference Manual.
 *
 * If the difference of the pressure potential between two cells is below the threshold
 * pressure, the pressure potential difference is assumed to be zero, if it is larger
 * than the threshold pressure, it is reduced by the threshold pressure.
 */
template <class TypeTag>
class EclThresholdPressure : public EclGenericThresholdPressure<GetPropType<TypeTag, Properties::Grid>,
                                                                GetPropType<TypeTag, Properties::GridView>,
                                                                GetPropType<TypeTag, Properties::ElementMapper>,
                                                                GetPropType<TypeTag, Properties::Scalar>>
{
    using BaseType = EclGenericThresholdPressure<GetPropType<TypeTag, Properties::Grid>,
                                                 GetPropType<TypeTag, Properties::GridView>,
                                                 GetPropType<TypeTag, Properties::ElementMapper>,
                                                 GetPropType<TypeTag, Properties::Scalar>>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ExtensiveQuantities = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    enum { dimWorld = GridView::dimensionworld };

    enum { enableExperiments = getPropValue<TypeTag, Properties::EnableExperiments>() };
    enum { numPhases = FluidSystem::numPhases };

public:
    EclThresholdPressure(const Simulator& simulator)
        : BaseType(simulator.vanguard().cartesianIndexMapper(),
                   simulator.vanguard().gridView(),
                   simulator.model().elementMapper(),
                   simulator.vanguard().eclState(),
                   simulator.vanguard().deck(),
                   enableExperiments)
        , simulator_(simulator)
    {
    }

    /*!
     * \brief Actually compute the threshold pressures over a face as a pre-compute step.
     */
    void finishInit()
    {
        this->BaseType::finishInit();
        if (this->enableThresholdPressure_ && !this->thpresDefault_.empty()) {
            this->computeDefaultThresholdPressures_();
            this->applyExplicitThresholdPressures_();
        }
    }

private:
    template<class Face,class Stencil,class ElemCtx>
    double calculateMaxDp(Face& face, Stencil& stencil,
                          ElemCtx& elemCtx,const unsigned& scvfIdx,
                          const unsigned& i,const unsigned& j,const unsigned& insideElemIdx,const unsigned& outsideElemIdx){
        typedef MathToolbox<Evaluation> Toolbox;
        elemCtx.updateIntensiveQuantities(/*timeIdx=*/0);
        elemCtx.updateExtensiveQuantities(/*timeIdx=*/0);
        // determine the maximum difference of the pressure of any phase over the
        // intersection
        Scalar pth = 0.0;
        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, /*timeIdx=*/0);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            unsigned upIdx = extQuants.upstreamIndex(phaseIdx);
            const auto& up = elemCtx.intensiveQuantities(upIdx, /*timeIdx=*/0);
            
            if (up.mobility(phaseIdx) > 0.0) {
                Scalar phaseVal = Toolbox::value(extQuants.pressureDifference(phaseIdx));
                pth = std::max(pth, std::abs(phaseVal));
            }
        }
        return pth;
    }
    
    template<class Face,class Stencil>
    double calculateMaxDp(Face& face, Stencil& stencil,
                          SmallElementContext<TypeTag>& elemCtx,const unsigned& scvfIdx,
                          const unsigned& i,const unsigned& j,
                          const unsigned& insideElemIdx,const unsigned& outsideElemIdx){        
        typedef MathToolbox<Evaluation> Toolbox;
        // determine the maximum difference of the pressure of any phase over the
        // intersection
        Scalar pth = 0.0;
        //const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, /*timeIdx=*/0);
        
        Scalar Vin = elemCtx.dofVolume(i, /*timeIdx=*/0);
        Scalar Vex = elemCtx.dofVolume(j, /*timeIdx=*/0);
        
        Scalar thpres = 0.0;//NB ??problem.thresholdPressure(globalIndexIn, globalIndexEx);
        
        // estimate the gravity correction: for performance reasons we use a simplified
        // approach for this flux module that assumes that gravity is constant and always
        // acts into the downwards direction. (i.e., no centrifuge experiments, sorry.)
        const auto& problem = elemCtx.problem();
        Scalar g = problem.gravity()[dimWorld - 1];
        
        const auto& intQuantsIn = elemCtx.intensiveQuantities(i, /*timeIdx*/0);
        const auto& intQuantsEx = elemCtx.intensiveQuantities(j, /*timeIdx*/0);
        
        // this is quite hacky because the dune grid interface does not provide a
        // cellCenterDepth() method (so we ask the problem to provide it). The "good"
        // solution would be to take the Z coordinate of the element centroids, but since
        // ECL seems to like to be inconsistent on that front, it needs to be done like
        // here...
        Scalar zIn = problem.dofCenterDepth(elemCtx, i, /*timeIdx*/0);
        Scalar zEx = problem.dofCenterDepth(elemCtx, j, /*timeIdx*/0);
        
        // the distances from the DOF's depths. (i.e., the additional depth of the
        // exterior DOF)
        Scalar distZ = zIn - zEx;
        
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            short dnIdx;
            //
            short upIdx;
            Evaluation pressureDifference;
            ExtensiveQuantities::calculatePhasePressureDiff_(upIdx,
                                                             dnIdx,
                                                             pressureDifference,
                                                             intQuantsIn,
                                                             intQuantsEx,
                                                             scvfIdx,//input
                                                             /*timeIdx*/0,//input
                                                             phaseIdx,//input
                                                             i,//input
                                                             j,//intput
                                                             Vin,
                                                             Vex,
                                                             insideElemIdx,
                                                             outsideElemIdx,
                                                             distZ*g,
                                                             thpres);
            const IntensiveQuantities& up = (upIdx == i) ? intQuantsIn : intQuantsEx;
            if (up.mobility(phaseIdx) > 0.0) {
                Scalar phaseVal = Toolbox::value(pressureDifference);
                pth = std::max(pth, std::abs(phaseVal));
            }
        }
        return pth;
    }
    // compute the defaults of the threshold pressures using the initial condition
    void computeDefaultThresholdPressures_()
    {
        const auto& vanguard = simulator_.vanguard();
        const auto& gridView = vanguard.gridView();

        typedef MathToolbox<Evaluation> Toolbox;
        // loop over the whole grid and compute the maximum gravity adjusted pressure
        // difference between two EQUIL regions.
        auto elemIt = gridView.template begin</*codim=*/ 0>();
        const auto& elemEndIt = gridView.template end</*codim=*/ 0>();
        ElementContext elemCtx(simulator_);
        simulator_.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
        for (; elemIt != elemEndIt; ++elemIt) {

            const auto& elem = *elemIt;
            if (elem.partitionType() != Dune::InteriorEntity)
                continue;

            elemCtx.updateStencil(elem);
            //
            
            const auto& stencil = elemCtx.stencil(/*timeIdx=*/0);
            

            for (unsigned scvfIdx = 0; scvfIdx < stencil.numInteriorFaces(); ++ scvfIdx) {
                
                const auto& face = stencil.interiorFace(scvfIdx);

                unsigned i = face.interiorIndex();
                unsigned j = face.exteriorIndex();

                unsigned insideElemIdx = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
                unsigned outsideElemIdx = elemCtx.globalSpaceIndex(j, /*timeIdx=*/0);

                unsigned equilRegionInside = this->elemEquilRegion_[insideElemIdx];
                unsigned equilRegionOutside = this->elemEquilRegion_[outsideElemIdx];
                
                if (equilRegionInside == equilRegionOutside)
                    // the current face is not at the boundary between EQUIL regions!
                    continue;
                const auto& problem = elemCtx.problem();
                // don't include connections with negligible flow
                const Evaluation& trans = problem.transmissibility(elemCtx, i, j);
                Scalar faceArea = face.area();
                if (std::abs(faceArea*getValue(trans)) < 1e-18)
                    continue;
                
                double pth = calculateMaxDp(face, stencil, elemCtx, scvfIdx,
                                            i, j,
                                            insideElemIdx, outsideElemIdx);
                // don't include connections with negligible flow

                int offset1 = equilRegionInside*this->numEquilRegions_ + equilRegionOutside;
                int offset2 = equilRegionOutside*this->numEquilRegions_ + equilRegionInside;

                this->thpresDefault_[offset1] = std::max(this->thpresDefault_[offset1], pth);
                this->thpresDefault_[offset2] = std::max(this->thpresDefault_[offset2], pth);
            }
        }

        // make sure that the threshold pressures is consistent for parallel
        // runs. (i.e. take the maximum of all processes)
        for (unsigned i = 0; i < this->thpresDefault_.size(); ++i)
            this->thpresDefault_[i] = gridView.comm().max(this->thpresDefault_[i]);
    }

    const Simulator& simulator_;
};

} // namespace Opm

#endif
