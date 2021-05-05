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

    void initParamsForElements(const EclipseState& eclState, size_t numElems)
    {
        const auto& fp = eclState.fieldProps();
        const auto& tableManager = eclState.getTableManager();
        bool has_heatcr = fp.has_double("HEATCR");
        bool has_thconr = fp.has_double("THCONR");
        bool has_thc = fp.has_double("THCROCK") || fp.has_double("THCOIL") || fp.has_double("THCGAS") || fp.has_double("THCWATER");

        if (has_heatcr)
            initHeatcr_(eclState, numElems);
        else if (tableManager.hasTables("SPECROCK"))
            initSpecrock_(eclState, numElems);
        else
            initNullRockEnergy_();

        if (has_thconr)
            initThconr_(eclState, numElems);
        else if (has_thc)
            initThc_(eclState, numElems);
        else
            initNullCond_();
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
    void initHeatcr_(const EclipseState& eclState,
                     size_t numElems)
    {
        solidEnergyApproach_ = SolidEnergyLawParams::heatcrApproach;
        // actually the value of the reference temperature does not matter for energy
        // conservation. We set it anyway to faciliate comparisons with ECL
        HeatcrLawParams::setReferenceTemperature(FluidSystem::surfaceTemperature);

        const auto& fp = eclState.fieldProps();
        const std::vector<double>& heatcrData  = fp.get_double("HEATCR");
        const std::vector<double>& heatcrtData = fp.get_double("HEATCRT");
        solidEnergyLawParams_.resize(numElems);
        for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            auto& elemParam = solidEnergyLawParams_[elemIdx];
            elemParam.setSolidEnergyApproach(SolidEnergyLawParams::heatcrApproach);
            auto& heatcrElemParams = elemParam.template getRealParams<SolidEnergyLawParams::heatcrApproach>();

            heatcrElemParams.setReferenceRockHeatCapacity(heatcrData[elemIdx]);
            heatcrElemParams.setDRockHeatCapacity_dT(heatcrtData[elemIdx]);
            heatcrElemParams.finalize();
            elemParam.finalize();
        }
    }

    /*!
     * \brief Initialize the parameters for the solid energy law using using SPECROCK and friends.
     */
    void initSpecrock_(const EclipseState& eclState,
                       size_t numElems)
    {
        solidEnergyApproach_ = SolidEnergyLawParams::specrockApproach;

        // initialize the element index -> SATNUM index mapping
        const auto& fp = eclState.fieldProps();
        const std::vector<int>& satnumData = fp.get_int("SATNUM");
        elemToSatnumIdx_.resize(numElems);
        for (unsigned elemIdx = 0; elemIdx < numElems; ++ elemIdx) {
            // satnumData contains Fortran-style indices, i.e., they start with 1 instead
            // of 0!
            elemToSatnumIdx_[elemIdx] = satnumData[elemIdx] - 1;
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
    void initNullRockEnergy_()
    {
        solidEnergyApproach_ = SolidEnergyLawParams::nullApproach;

        solidEnergyLawParams_.resize(1);
        solidEnergyLawParams_[0].finalize();
    }

    /*!
     * \brief Initialize the parameters for the thermal conduction law using THCONR and friends.
     */
    void initThconr_(const EclipseState& eclState,
                     size_t numElems)
    {
        thermalConductivityApproach_ = ThermalConductionLawParams::thconrApproach;

        const auto& fp = eclState.fieldProps();
        std::vector<double> thconrData;
        std::vector<double> thconsfData;
        if (fp.has_double("THCONR"))
            thconrData  = fp.get_double("THCONR");

        if (fp.has_double("THCONSF"))
            thconsfData = fp.get_double("THCONSF");

        thermalConductionLawParams_.resize(numElems);
        for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            auto& elemParams = thermalConductionLawParams_[elemIdx];
            elemParams.setThermalConductionApproach(ThermalConductionLawParams::thconrApproach);
            auto& thconrElemParams = elemParams.template getRealParams<ThermalConductionLawParams::thconrApproach>();

            double thconr = thconrData.empty()   ? 0.0 : thconrData[elemIdx];
            double thconsf = thconsfData.empty() ? 0.0 : thconsfData[elemIdx];
            thconrElemParams.setReferenceTotalThermalConductivity(thconr);
            thconrElemParams.setDTotalThermalConductivity_dSg(thconsf);

            thconrElemParams.finalize();
            elemParams.finalize();
        }
    }

    /*!
     * \brief Initialize the parameters for the thermal conduction law using THCROCK and friends.
     */
    void initThc_(const EclipseState& eclState,
                  size_t numElems)
    {
        thermalConductivityApproach_ = ThermalConductionLawParams::thcApproach;

        const auto& fp = eclState.fieldProps();
        std::vector<double> thcrockData;
        std::vector<double> thcoilData;
        std::vector<double> thcgasData;
        std::vector<double> thcwaterData = fp.get_double("THCWATER");

        if (fp.has_double("THCROCK"))
            thcrockData = fp.get_double("THCROCK");

        if (fp.has_double("THCOIL"))
            thcoilData = fp.get_double("THCOIL");

        if (fp.has_double("THCGAS"))
            thcgasData = fp.get_double("THCGAS");

        if (fp.has_double("THCWATER"))
            thcwaterData = fp.get_double("THCWATER");

        const std::vector<double>& poroData = fp.get_double("PORO");

        thermalConductionLawParams_.resize(numElems);
        for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            auto& elemParams = thermalConductionLawParams_[elemIdx];
            elemParams.setThermalConductionApproach(ThermalConductionLawParams::thcApproach);
            auto& thcElemParams = elemParams.template getRealParams<ThermalConductionLawParams::thcApproach>();

            thcElemParams.setPorosity(poroData[elemIdx]);
            double thcrock = thcrockData.empty()    ? 0.0 : thcrockData[elemIdx];
            double thcoil = thcoilData.empty()      ? 0.0 : thcoilData[elemIdx];
            double thcgas = thcgasData.empty()      ? 0.0 : thcgasData[elemIdx];
            double thcwater = thcwaterData.empty()  ? 0.0 : thcwaterData[elemIdx];
            thcElemParams.setThcrock(thcrock);
            thcElemParams.setThcoil(thcoil);
            thcElemParams.setThcgas(thcgas);
            thcElemParams.setThcwater(thcwater);

            thcElemParams.finalize();
            elemParams.finalize();
        }
    }

    /*!
     * \brief Disable thermal conductivity
     */
    void initNullCond_()
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
