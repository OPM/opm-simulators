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
 * \copydoc Opm::EclThermalLawManager
 */
#if ! HAVE_ECL_INPUT
#error "Eclipse input support in opm-common is required to use the ECL thermal law manager!"
#endif

#ifndef OPM_ECL_THERMAL_LAW_MANAGER_HPP
#define OPM_ECL_THERMAL_LAW_MANAGER_HPP

#include "EclSolidEnergyLawMultiplexer.hpp"
#include "EclSolidEnergyLawMultiplexerParams.hpp"

#include "EclThermalConductionLawMultiplexer.hpp"
#include "EclThermalConductionLawMultiplexerParams.hpp"

#include <opm/material/common/Exceptions.hpp>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

namespace Opm {

/*!
 * \ingroup fluidmatrixinteractions
 *
 * \brief Provides an simple way to create and manage the thermal law objects
 *        for a complete ECL deck.
 */
template <class Scalar, class FluidSystem>
class EclThermalLawManager
{
public:
    typedef EclSolidEnergyLawMultiplexer<Scalar, FluidSystem> SolidEnergyLaw;
    typedef typename SolidEnergyLaw::Params SolidEnergyLawParams;
    typedef typename SolidEnergyLawParams::HeatcrLawParams HeatcrLawParams;
    typedef typename SolidEnergyLawParams::SpecrockLawParams SpecrockLawParams;

    typedef EclThermalConductionLawMultiplexer<Scalar, FluidSystem> ThermalConductionLaw;
    typedef typename ThermalConductionLaw::Params ThermalConductionLawParams;

    EclThermalLawManager()
    {
        solidEnergyApproach_ = SolidEnergyLawParams::undefinedApproach;
        thermalConductivityApproach_ = ThermalConductionLawParams::undefinedApproach;
    }

    void initFromDeck(const Opm::Deck& deck,
                      const Opm::EclipseState& eclState,
                      const std::vector<int>& compressedToCartesianElemIdx)
    {
        const auto& props = eclState.get3DProperties();
        if (props.hasDeckDoubleGridProperty("HEATCR"))
            initHeatcr_(deck, eclState, compressedToCartesianElemIdx);
        else if (deck.hasKeyword("SPECROCK"))
            initSpecrock_(deck, eclState, compressedToCartesianElemIdx);
        else
            initNullRockEnergy_(deck, eclState, compressedToCartesianElemIdx);

        if (props.hasDeckDoubleGridProperty("THCONR"))
            initThconr_(deck, eclState, compressedToCartesianElemIdx);
        else if (props.hasDeckDoubleGridProperty("THCROCK")
                 || props.hasDeckDoubleGridProperty("THCOIL")
                 || props.hasDeckDoubleGridProperty("THCGAS")
                 || props.hasDeckDoubleGridProperty("THCWATER"))
            initThc_(deck, eclState, compressedToCartesianElemIdx);
        else
            initNullCond_(deck, eclState, compressedToCartesianElemIdx);
    }

    const SolidEnergyLawParams& solidEnergyLawParams(unsigned elemIdx) const
    {
        switch (solidEnergyApproach_) {
        case SolidEnergyLawParams::heatcrApproach:
            assert(0 <= elemIdx && elemIdx <  solidEnergyLawParams_.size());
            return solidEnergyLawParams_[elemIdx];

        case SolidEnergyLawParams::specrockApproach:
        {
            assert(0 <= elemIdx && elemIdx <  elemToSatnumIdx_.size());
            unsigned satnumIdx = elemToSatnumIdx_[elemIdx];
            assert(0 <= satnumIdx && satnumIdx <  solidEnergyLawParams_.size());
            return solidEnergyLawParams_[satnumIdx];
        }

        case SolidEnergyLawParams::nullApproach:
            return solidEnergyLawParams_[0];

        default:
            throw std::runtime_error("Attempting to retrieve solid energy storage parameters "
                                     "without a known approach being defined by the deck.");
        }
    }

    const ThermalConductionLawParams& thermalConductionLawParams(unsigned elemIdx) const
    {
        switch (thermalConductivityApproach_) {
        case ThermalConductionLawParams::thconrApproach:
        case ThermalConductionLawParams::thcApproach:
            assert(0 <= elemIdx && elemIdx <  thermalConductionLawParams_.size());
            return thermalConductionLawParams_[elemIdx];

        case ThermalConductionLawParams::nullApproach:
            return thermalConductionLawParams_[0];

        default:
            throw std::runtime_error("Attempting to retrieve thermal conduction parameters without "
                                     "a known approach being defined by the deck.");
        }
    }

private:
    /*!
     * \brief Initialize the parameters for the solid energy law using using HEATCR and friends.
     */
    void initHeatcr_(const Opm::Deck& deck,
                     const Opm::EclipseState& eclState,
                     const std::vector<int>& compressedToCartesianElemIdx)
    {
        solidEnergyApproach_ = SolidEnergyLawParams::heatcrApproach;

        const auto& props = eclState.get3DProperties();

        const std::vector<double>& heatcrData = props.getDoubleGridProperty("HEATCR").getData();
        const std::vector<double>& heatcrtData = props.getDoubleGridProperty("HEATCRT").getData();

        // actually the value of the reference temperature does not matter for energy
        // conservation. We set it anyway to faciliate comparisons with ECL
        HeatcrLawParams::setReferenceTemperature(FluidSystem::surfaceTemperature);

        unsigned numElems = compressedToCartesianElemIdx.size();
        solidEnergyLawParams_.resize(numElems);
        for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            int cartElemIdx = compressedToCartesianElemIdx[elemIdx];

            auto& elemParam = solidEnergyLawParams_[elemIdx];

            elemParam.setSolidEnergyApproach(SolidEnergyLawParams::heatcrApproach);

            auto& heatcrElemParams = elemParam.template getRealParams<SolidEnergyLawParams::heatcrApproach>();
            heatcrElemParams.setReferenceRockHeatCapacity(heatcrData[cartElemIdx]);
            heatcrElemParams.setDRockHeatCapacity_dT(heatcrtData[cartElemIdx]);
            heatcrElemParams.finalize();

            elemParam.finalize();
        }
    }

    /*!
     * \brief Initialize the parameters for the solid energy law using using SPECROCK and friends.
     */
    void initSpecrock_(const Opm::Deck& deck,
                       const Opm::EclipseState& eclState,
                       const std::vector<int>& compressedToCartesianElemIdx)
    {
        solidEnergyApproach_ = SolidEnergyLawParams::specrockApproach;

        // initialize the element index -> SATNUM index mapping
        const auto& props = eclState.get3DProperties();
        const std::vector<int>& satnumData = props.getIntGridProperty("SATNUM").getData();
        elemToSatnumIdx_.resize(compressedToCartesianElemIdx.size());
        for (unsigned elemIdx = 0; elemIdx < compressedToCartesianElemIdx.size(); ++ elemIdx) {
            unsigned cartesianElemIdx = compressedToCartesianElemIdx[elemIdx];

            // satnumData contains Fortran-style indices, i.e., they start with 1 instead
            // of 0!
            elemToSatnumIdx_[elemIdx] = satnumData[cartesianElemIdx] - 1;
        }

        // internalize the SPECROCK table
        unsigned numSatRegions = eclState.runspec().tabdims().getNumSatTables();
        const auto& tableManager = eclState.getTableManager();
        solidEnergyLawParams_.resize(numSatRegions);
        for (unsigned satnumIdx = 0; satnumIdx < numSatRegions; ++satnumIdx) {
            const auto& specrockTable = tableManager.getSpecrockTables()[satnumIdx];

            auto& multiplexerParams = solidEnergyLawParams_[satnumIdx];

            multiplexerParams.setSolidEnergyApproach(SolidEnergyLawParams::specrockApproach);

            auto& specrockParams = multiplexerParams.template getRealParams<SolidEnergyLawParams::specrockApproach>();
            const auto& temperatureColumn = specrockTable.getColumn("TEMPERATURE");
            const auto& cvRockColumn = specrockTable.getColumn("CV_ROCK");
            specrockParams.setHeatCapacities(temperatureColumn, cvRockColumn);
            specrockParams.finalize();

            multiplexerParams.finalize();
        }
    }

    /*!
     * \brief Specify the solid energy law by setting heat capacity of rock to 0
     */
    void initNullRockEnergy_(const Opm::Deck& deck,
                      const Opm::EclipseState& eclState,
                      const std::vector<int>& compressedToCartesianElemIdx)
    {
        solidEnergyApproach_ = SolidEnergyLawParams::nullApproach;

        solidEnergyLawParams_.resize(1);
        solidEnergyLawParams_[0].finalize();
    }

    /*!
     * \brief Initialize the parameters for the thermal conduction law using THCONR and friends.
     */
    void initThconr_(const Opm::Deck& deck,
                     const Opm::EclipseState& eclState,
                     const std::vector<int>& compressedToCartesianElemIdx)
    {
        thermalConductivityApproach_ = ThermalConductionLawParams::thconrApproach;

        const auto& props = eclState.get3DProperties();

        const std::vector<double>& thconrData = props.getDoubleGridProperty("THCONR").getData();
        const std::vector<double>& thconsfData = props.getDoubleGridProperty("THCONSF").getData();

        unsigned numElems = compressedToCartesianElemIdx.size();
        thermalConductionLawParams_.resize(numElems);
        for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            int cartElemIdx = compressedToCartesianElemIdx[elemIdx];

            auto& elemParams = thermalConductionLawParams_[elemIdx];

            elemParams.setThermalConductionApproach(ThermalConductionLawParams::thconrApproach);

            auto& thconrElemParams = elemParams.template getRealParams<ThermalConductionLawParams::thconrApproach>();
            thconrElemParams.setReferenceTotalThermalConductivity(thconrData[cartElemIdx]);
            thconrElemParams.setDTotalThermalConductivity_dSg(thconsfData[cartElemIdx]);
            thconrElemParams.finalize();

            elemParams.finalize();
        }
    }

    /*!
     * \brief Initialize the parameters for the thermal conduction law using THCROCK and friends.
     */
    void initThc_(const Opm::Deck& deck,
                  const Opm::EclipseState& eclState,
                  const std::vector<int>& compressedToCartesianElemIdx)
    {
        thermalConductivityApproach_ = ThermalConductionLawParams::thcApproach;

        const auto& props = eclState.get3DProperties();

        const std::vector<double>& thcrockData = props.getDoubleGridProperty("THCROCK").getData();
        const std::vector<double>& thcoilData = props.getDoubleGridProperty("THCOIL").getData();
        const std::vector<double>& thcgasData = props.getDoubleGridProperty("THCGAS").getData();
        const std::vector<double>& thcwaterData = props.getDoubleGridProperty("THCWATER").getData();
        const std::vector<double>& poroData = props.getDoubleGridProperty("PORO").getData();

        unsigned numElems = compressedToCartesianElemIdx.size();
        thermalConductionLawParams_.resize(numElems);
        for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            int cartElemIdx = compressedToCartesianElemIdx[elemIdx];

            auto& elemParams = thermalConductionLawParams_[elemIdx];

            elemParams.setThermalConductionApproach(ThermalConductionLawParams::thcApproach);

            auto& thcElemParams = elemParams.template getRealParams<ThermalConductionLawParams::thcApproach>();
            thcElemParams.setPorosity(poroData[cartElemIdx]);
            thcElemParams.setThcrock(thcrockData[cartElemIdx]);
            thcElemParams.setThcoil(thcoilData[cartElemIdx]);
            thcElemParams.setThcgas(thcgasData[cartElemIdx]);
            thcElemParams.setThcwater(thcwaterData[cartElemIdx]);
            thcElemParams.finalize();

            elemParams.finalize();
        }
    }

    /*!
     * \brief Disable thermal conductivity
     */
    void initNullCond_(const Opm::Deck& deck,
                       const Opm::EclipseState& eclState,
                       const std::vector<int>& compressedToCartesianElemIdx)
    {
        thermalConductivityApproach_ = ThermalConductionLawParams::nullApproach;

        thermalConductionLawParams_.resize(1);
        thermalConductionLawParams_[0].finalize();
    }

private:
    typename ThermalConductionLawParams::ThermalConductionApproach thermalConductivityApproach_;
    typename SolidEnergyLawParams::SolidEnergyApproach solidEnergyApproach_;

    std::vector<unsigned> elemToSatnumIdx_;

    std::vector<SolidEnergyLawParams> solidEnergyLawParams_;
    std::vector<ThermalConductionLawParams> thermalConductionLawParams_;
};
} // namespace Opm

#endif
