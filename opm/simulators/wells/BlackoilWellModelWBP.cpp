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

#include <config.h>
#include <opm/simulators/wells/BlackoilWellModelWBP.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>

#include <cassert>

namespace Opm {

template<typename Scalar, typename IndexTraits>
BlackoilWellModelWBP<Scalar, IndexTraits>::
BlackoilWellModelWBP(BlackoilWellModelGeneric<Scalar, IndexTraits>& well_model)
    : well_model_(well_model)
    , wbpCalculationService_(well_model.eclipseState().gridDims(), well_model.comm())
{}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelWBP<Scalar, IndexTraits>::
initializeSources(typename ParallelWBPCalculation<Scalar>::GlobalToLocal index,
                  typename ParallelWBPCalculation<Scalar>::Evaluator eval)
{
    this->wbpCalculationService_.localCellIndex(index).evalCellSource(eval);
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelWBP<Scalar, IndexTraits>::
registerOpenWellsForWBPCalculation()
{
    assert(this->wbpCalcMap_.size() ==
           static_cast<std::size_t>(well_model_.numLocalWells()));

    for (auto& wbpCalc : this->wbpCalcMap_) {
        wbpCalc.openWellIdx_.reset();
    }

    auto openWellIdx = typename std::vector<WellInterfaceGeneric<Scalar, IndexTraits>*>::size_type{0};
    for (const auto* openWell : well_model_.genericWells()) {
        this->wbpCalcMap_[openWell->indexOfWell()].openWellIdx_ = openWellIdx++;
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelWBP<Scalar, IndexTraits>::
initializeWBPCalculationService()
{
    this->wbpCalcMap_.clear();
    this->wbpCalcMap_.resize(well_model_.numLocalWells());

    this->registerOpenWellsForWBPCalculation();

    auto wellID = std::size_t{0};
    for (const auto& well : well_model_.eclWells()) {
        this->wbpCalcMap_[wellID].wbpCalcIdx_ = this->wbpCalculationService_
            .createCalculator(well,
                              well_model_.parallelWellInfo(wellID),
                              well_model_.connectionIndexMap(wellID).local(),
                              this->makeWellSourceEvaluatorFactory(wellID));

        ++wellID;
    }

    this->wbpCalculationService_.defineCommunication();
}

template<typename Scalar, typename IndexTraits>
data::WellBlockAveragePressures
BlackoilWellModelWBP<Scalar, IndexTraits>::
computeWellBlockAveragePressures(const Scalar gravity) const
{
    auto wbpResult = data::WellBlockAveragePressures{};

    using Calculated = typename PAvgCalculatorResult<Scalar>::WBPMode;
    using Output = data::WellBlockAvgPress::Quantity;

    this->wbpCalculationService_.collectDynamicValues();

    const auto numWells = well_model_.numLocalWells();
    for (auto wellID = 0*numWells; wellID < numWells; ++wellID) {
        const auto calcIdx = this->wbpCalcMap_[wellID].wbpCalcIdx_;
        const auto& well = well_model_.eclWells()[wellID];

        if (! well.hasRefDepth()) {
            // Can't perform depth correction without at least a
            // fall-back datum depth.
            continue;
        }

        this->wbpCalculationService_
            .inferBlockAveragePressures(calcIdx, well.pavg(),
                                        gravity,
                                        well.getWPaveRefDepth());

        const auto& result = this->wbpCalculationService_
            .averagePressures(calcIdx);

        auto& reported = wbpResult.values[well.name()];

        reported[Output::WBP]  = result.value(Calculated::WBP);
        reported[Output::WBP4] = result.value(Calculated::WBP4);
        reported[Output::WBP5] = result.value(Calculated::WBP5);
        reported[Output::WBP9] = result.value(Calculated::WBP9);
    }

    return wbpResult;
}

template<typename Scalar, typename IndexTraits>
typename ParallelWBPCalculation<Scalar>::EvaluatorFactory
BlackoilWellModelWBP<Scalar, IndexTraits>::
makeWellSourceEvaluatorFactory(const std::vector<Well>::size_type wellIdx) const
{
    using Span = typename PAvgDynamicSourceData<Scalar>::template SourceDataSpan<Scalar>;
    using Item = typename Span::Item;

    return [wellIdx, this]() -> typename ParallelWBPCalculation<Scalar>::Evaluator
    {
        if (! this->wbpCalcMap_[wellIdx].openWellIdx_.has_value()) {
            // Well is stopped/shut.  Return evaluator for stopped wells.
            return []([[maybe_unused]] const int connIdx, Span sourceTerm)
            {
                // Well/connection is stopped/shut.  Set all items to
                // zero.

                sourceTerm
                    .set(Item::Pressure      , 0.0)
                    .set(Item::PoreVol       , 0.0)
                    .set(Item::MixtureDensity, 0.0)
                    .set(Item::Depth         , 0.0)
                    ;
            };
        }

        // Well is open.  Return an evaluator for open wells/open connections.
        return [this, wellPtr = well_model_.genericWells()[*this->wbpCalcMap_[wellIdx].openWellIdx_]]
            (const int connIdx, Span sourceTerm)
        {
            // Note: The only item which actually matters for the WBP
            // calculation at the well reservoir connection level is the
            // mixture density.  Set other items to zero.

            const auto& connIdxMap =
                well_model_.connectionIndexMap(wellPtr->indexOfWell());

            const auto rho = wellPtr->
                connectionDensity(connIdxMap.global(connIdx),
                                  connIdxMap.open(connIdx));

            sourceTerm
                .set(Item::Pressure      , 0.0)
                .set(Item::PoreVol       , 0.0)
                .set(Item::MixtureDensity, rho)
                .set(Item::Depth         , 0.0)
                ;
        };
    };
}

template class BlackoilWellModelWBP<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class BlackoilWellModelWBP<float, BlackOilDefaultFluidSystemIndices>;
#endif

}
