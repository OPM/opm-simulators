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
 * \copydoc Opm::EclEpsTwoPhaseLawPoints
 */
#ifndef OPM_ECL_EPS_SCALING_POINTS_HPP
#define OPM_ECL_EPS_SCALING_POINTS_HPP

#include "EclEpsConfig.hpp"

#if HAVE_ECL_INPUT
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Deck/DeckRecord.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/GridProperty.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SgfnTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SgofTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SlgofTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Sof2Table.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Sof3Table.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SwfnTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SwofTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#endif

#include <opm/material/common/Means.hpp>

#include <array>
#include <string>
#include <iostream>
#include <cassert>
#include <algorithm>

namespace Opm {
/*!
 * \brief Collects all grid properties which are relevant for end point scaling.
 *
 * This class is used for both, the drainage and the imbibition variants of the ECL
 * keywords.
 */
class EclEpsGridProperties
{
    typedef std::vector<int> IntData;
    typedef std::vector<double> DoubleData;

public:
#if HAVE_ECL_INPUT
    void initFromDeck(const Opm::Deck& /* deck */,
                      const Opm::EclipseState& eclState,
                      bool useImbibition)
    {
        std::string kwPrefix = useImbibition?"I":"";

        if (useImbibition)
            satnum = &eclState.get3DProperties().getIntGridProperty("IMBNUM").getData();
        else
            satnum = &eclState.get3DProperties().getIntGridProperty("SATNUM").getData();

        retrieveGridPropertyData_(&swl, eclState, kwPrefix+"SWL");
        retrieveGridPropertyData_(&sgl, eclState, kwPrefix+"SGL");
        retrieveGridPropertyData_(&swcr, eclState, kwPrefix+"SWCR");
        retrieveGridPropertyData_(&sgcr, eclState, kwPrefix+"SGCR");
        retrieveGridPropertyData_(&sowcr, eclState, kwPrefix+"SOWCR");
        retrieveGridPropertyData_(&sogcr, eclState, kwPrefix+"SOGCR");
        retrieveGridPropertyData_(&swu, eclState, kwPrefix+"SWU");
        retrieveGridPropertyData_(&sgu, eclState, kwPrefix+"SGU");
        retrieveGridPropertyData_(&pcw, eclState, kwPrefix+"PCW");
        retrieveGridPropertyData_(&pcg, eclState, kwPrefix+"PCG");
        retrieveGridPropertyData_(&krw, eclState, kwPrefix+"KRW");
        retrieveGridPropertyData_(&kro, eclState, kwPrefix+"KRO");
        retrieveGridPropertyData_(&krg, eclState, kwPrefix+"KRG");

        // _may_ be needed to calculate the Leverett capillary pressure scaling factor
        const auto& ecl3dProps = eclState.get3DProperties();
        poro = &ecl3dProps.getDoubleGridProperty("PORO").getData();

        if (ecl3dProps.hasDeckDoubleGridProperty("PERMX")) {
            permx = &ecl3dProps.getDoubleGridProperty("PERMX").getData();
            permy = permx;
            permz = permx;
        }

        if (ecl3dProps.hasDeckDoubleGridProperty("PERMY"))
            permy = &ecl3dProps.getDoubleGridProperty("PERMY").getData();

        if (ecl3dProps.hasDeckDoubleGridProperty("PERMZ"))
            permz = &ecl3dProps.getDoubleGridProperty("PERMZ").getData();
    }
#endif

    const IntData* satnum;

    const DoubleData* swl;
    const DoubleData* sgl;
    const DoubleData* swcr;
    const DoubleData* sgcr;
    const DoubleData* sowcr;
    const DoubleData* sogcr;
    const DoubleData* swu;
    const DoubleData* sgu;
    const DoubleData* pcw;
    const DoubleData* pcg;
    const DoubleData* krw;
    const DoubleData* kro;
    const DoubleData* krg;
    const DoubleData* poro;
    const DoubleData* permx;
    const DoubleData* permy;
    const DoubleData* permz;

private:
#if HAVE_ECL_INPUT
    // this method makes sure that a grid property is not created if it is not explicitly
    // mentioned in the deck. (saves memory.)
    void retrieveGridPropertyData_(const DoubleData **data,
                                   const Opm::EclipseState& eclState,
                                   const std::string& properyName)
    {
        (*data) = 0;
        if (eclState.get3DProperties().hasDeckDoubleGridProperty(properyName))
            (*data) = &eclState.get3DProperties().getDoubleGridProperty(properyName).getData();
    }
#endif
};

/*!
 * \brief This structure represents all values which can be possibly used as scaling
 *        points in the endpoint scaling code.
 *
 * Depending on the exact configuration, some of these quantities are not used as actual
 * scaling points. It is easier to extract all of them at once, though.
 */
template <class Scalar>
struct EclEpsScalingPointsInfo
{
    // connate saturations
    Scalar Swl; // oil
    Scalar Sgl; // gas
    Scalar Sowl; // oil for the oil-water system
    Scalar Sogl; // oil for the gas-oil system

    // critical water and gas saturations
    Scalar krCriticalEps; // relative permeability below which a saturation is considered
                          // to be critical
    Scalar Swcr; // oil
    Scalar Sgcr; // gas
    Scalar Sowcr; // oil for the oil-water system
    Scalar Sogcr; // oil for the gas-oil system

    // maximum saturations
    Scalar Swu; // oil
    Scalar Sgu; // gas
    Scalar Sowu; // oil for the oil-water system
    Scalar Sogu; // oil for the gas-oil system

    // maximum capillary pressures
    Scalar maxPcow; // maximum capillary pressure of the oil-water system
    Scalar maxPcgo; // maximum capillary pressure of the gas-oil system

    // the Leverett capillary pressure scaling factors. (those only make sense for the
    // scaled points, for the unscaled ones they are 1.0.)
    Scalar pcowLeverettFactor;
    Scalar pcgoLeverettFactor;

    // maximum relative permabilities
    Scalar maxKrw; // maximum relative permability of water
    Scalar maxKrow; // maximum relative permability of oil in the oil-water system
    Scalar maxKrog; // maximum relative permability of oil in the gas-oil system
    Scalar maxKrg; // maximum relative permability of gas

    void print() const
    {
        std::cout << "    Swl: " << Swl << "\n"
                  << "    Sgl: " << Sgl << "\n"
                  << "    Sowl: " << Sowl << "\n"
                  << "    Sogl: " << Sogl << "\n"
                  << "    Swcr: " << Swcr << "\n"
                  << "    Sgcr: " << Sgcr << "\n"
                  << "    Sowcr: " << Sowcr << "\n"
                  << "    Sogcr: " << Sogcr << "\n"
                  << "    Swu: " << Swu << "\n"
                  << "    Sgu: " << Sgu << "\n"
                  << "    Sowu: " << Sowu << "\n"
                  << "    Sogu: " << Sogu << "\n"
                  << "    maxPcow: " << maxPcow << "\n"
                  << "    maxPcgo: " << maxPcgo << "\n"
                  << "    pcowLeverettFactor: " << pcowLeverettFactor << "\n"
                  << "    pcgoLeverettFactor: " << pcgoLeverettFactor << "\n"
                  << "    maxKrw: " << maxKrw << "\n"
                  << "    maxKrg: " << maxKrg << "\n"
                  << "    maxKrow: " << maxKrow << "\n"
                  << "    maxKrog: " << maxKrog << "\n";
    }

#if HAVE_ECL_INPUT
    /*!
     * \brief Extract the values of the unscaled scaling parameters.
     *
     * I.e., the values which are used for the nested Fluid-Matrix interactions and which
     * are produced by them.
     */
    void extractUnscaled(const Opm::Deck& deck,
                         const Opm::EclipseState& eclState,
                         unsigned satRegionIdx)
    {
        // determine the value of the relative permeability below which the corresponding
        // saturation is considered to be critical
        krCriticalEps = Opm::ParserKeywords::TOLCRIT::VALUE::defaultValue;
        if (deck.hasKeyword("TOLCRIT")) {
            const Opm::DeckKeyword& tolcritKeyword = deck.getKeyword("TOLCRIT");

            krCriticalEps = tolcritKeyword.getRecord(0).getItem("VALUE").getSIDouble(0);
        }

        const auto& tables = eclState.getTableManager();
        const TableContainer&  swofTables = tables.getSwofTables();
        const TableContainer&  sgofTables = tables.getSgofTables();
        const TableContainer& slgofTables = tables.getSlgofTables();
        const TableContainer&  swfnTables = tables.getSwfnTables();
        const TableContainer&  sgfnTables = tables.getSgfnTables();
        const TableContainer&  sof3Tables = tables.getSof3Tables();
        const TableContainer&  sof2Tables = tables.getSof2Tables();


        bool hasWater = deck.hasKeyword("WATER");
        bool hasGas = deck.hasKeyword("GAS");
        bool hasOil = deck.hasKeyword("OIL");

        if (!hasWater) {
            Swl = 0.0;
            Swu = 0.0;
            Swcr = 0.0;
            bool family1 = (!sgofTables.empty() || !slgofTables.empty());
            bool family2 = !sgfnTables.empty() && !sof2Tables.empty();
            if (family1) {
                if (!sgofTables.empty())
                    extractUnscaledSgof_(sgofTables.getTable<SgofTable>(satRegionIdx));
                else {
                    assert(!slgofTables.empty());
                    extractUnscaledSlgof_(slgofTables.getTable<SlgofTable>(satRegionIdx));
                }
            } else if (family2) {
                extractUnscaledSgfn_(sgfnTables.getTable<SgfnTable>(satRegionIdx));
                extractUnscaledSof2_(sof2Tables.getTable<Sof2Table>(satRegionIdx));
            }
            else {
                throw std::domain_error("No valid saturation keyword family specified");
            }
            return;
        }
        else if (!hasGas) {
            Sgl = 0.0;
            Sgu = 0.0;
            Sgcr = 0.0;
            bool family1 = !swofTables.empty();
            bool family2 = !swfnTables.empty() && !sof2Tables.empty();
            if (family1) {
                extractUnscaledSwof_(swofTables.getTable<SwofTable>(satRegionIdx));
            } else if (family2) {
                extractUnscaledSwfn_(swfnTables.getTable<SwfnTable>(satRegionIdx));
                extractUnscaledSof2_(sof2Tables.getTable<Sof2Table>(satRegionIdx));
            }
            else {
                throw std::domain_error("No valid saturation keyword family specified");
            }
            return;
        }

        bool family1 = (!sgofTables.empty() || !slgofTables.empty()) && !swofTables.empty();
        bool family2 = !swfnTables.empty() && !sgfnTables.empty() && !sof3Tables.empty();

        // so far, only water-oil and oil-gas simulations are supported, i.e.,
        // there's no gas-water yet.
        if (!hasWater || !hasGas || !hasOil)
            throw std::domain_error("The specified phase configuration is not suppored");

        if (family1) {
            extractUnscaledSwof_(swofTables.getTable<SwofTable>(satRegionIdx));

            if (!sgofTables.empty()) {
                // gas-oil parameters are specified using the SGOF keyword
                extractUnscaledSgof_(sgofTables.getTable<SgofTable>(satRegionIdx));
            }
            else {
                // gas-oil parameters are specified using the SLGOF keyword
                assert(!slgofTables.empty());

                extractUnscaledSlgof_(slgofTables.getTable<SlgofTable>(satRegionIdx));
            }
        }
        else if (family2) {
            extractUnscaledSwfn_(swfnTables.getTable<SwfnTable>(satRegionIdx));
            extractUnscaledSgfn_(sgfnTables.getTable<SgfnTable>(satRegionIdx));
            extractUnscaledSof3_(sof3Tables.getTable<Sof3Table>(satRegionIdx));
        }
        else {
            throw std::domain_error("No valid saturation keyword family specified");
        }

        // there are no "unscaled" Leverett factors, so we just set them to 1.0
        pcowLeverettFactor = 1.0;
        pcgoLeverettFactor = 1.0;
    }

    /*!
     * \brief Extract the values of the scaled scaling parameters.
     *
     * I.e., the values which are "seen" by the physical model.
     */
    void extractScaled(const Opm::EclipseState& eclState,
                       const EclEpsGridProperties& epsProperties,
                       unsigned cartesianCellIdx)
    {
        // overwrite the unscaled values with the values for the cell if it is
        // explicitly specified by the corresponding keyword.
        extractGridPropertyValue_(Swl, epsProperties.swl, cartesianCellIdx);
        extractGridPropertyValue_(Sgl, epsProperties.sgl, cartesianCellIdx);
        extractGridPropertyValue_(Swcr, epsProperties.swcr, cartesianCellIdx);
        extractGridPropertyValue_(Sgcr, epsProperties.sgcr, cartesianCellIdx);
        extractGridPropertyValue_(Sowcr, epsProperties.sowcr, cartesianCellIdx);
        extractGridPropertyValue_(Sogcr, epsProperties.sogcr, cartesianCellIdx);
        extractGridPropertyValue_(Swu, epsProperties.swu, cartesianCellIdx);
        extractGridPropertyValue_(Sgu, epsProperties.sgu, cartesianCellIdx);

        extractGridPropertyValue_(maxPcow, epsProperties.pcw, cartesianCellIdx);
        extractGridPropertyValue_(maxPcgo, epsProperties.pcg, cartesianCellIdx);

        extractGridPropertyValue_(maxKrw, epsProperties.krw, cartesianCellIdx);
        extractGridPropertyValue_(maxKrg, epsProperties.krg, cartesianCellIdx);

        // quite likely that's wrong!
        extractGridPropertyValue_(maxKrow, epsProperties.kro, cartesianCellIdx);
        extractGridPropertyValue_(maxKrog, epsProperties.kro, cartesianCellIdx);

        // compute the Leverett capillary pressure scaling factors if applicable.  note
        // that this needs to be done using non-SI units to make it correspond to the
        // documentation.
        pcowLeverettFactor = 1.0;
        pcgoLeverettFactor = 1.0;
        if (eclState.getTableManager().useJFunc()) {
            const auto& jfunc = eclState.getTableManager().getJFunc();
            const auto& jfuncDir = jfunc.direction();

            Scalar perm;
            if (jfuncDir == Opm::JFunc::Direction::X)
                perm =
                    (*epsProperties.permx)[cartesianCellIdx];
            else if (jfuncDir == Opm::JFunc::Direction::Y)
                perm =
                    (*epsProperties.permy)[cartesianCellIdx];
            else if (jfuncDir == Opm::JFunc::Direction::Z)
                perm =
                    (*epsProperties.permz)[cartesianCellIdx];
            else if (jfuncDir == Opm::JFunc::Direction::XY)
                // TODO: verify that this really is the arithmetic mean. (the
                // documentation just says that the "average" should be used, IMO the
                // harmonic mean would be more appropriate because that's what's usually
                // applied when calculating the fluxes.)
                perm =
                    Opm::arithmeticMean((*epsProperties.permx)[cartesianCellIdx],
                                        (*epsProperties.permy)[cartesianCellIdx]);
            else
                throw std::runtime_error("Illegal direction indicator for the JFUNC "
                                         "keyword ("+std::to_string(int(jfuncDir))+")");

            // convert permeability from m^2 to mD
            perm *= 1.01325e15;

            Scalar poro = (*epsProperties.poro)[cartesianCellIdx];
            Scalar alpha = jfunc.alphaFactor();
            Scalar beta = jfunc.betaFactor();

            // the part of the Leverett capillary pressure which does not depend on
            // surface tension.
            Scalar commonFactor = std::pow(poro, alpha)/std::pow(perm, beta);

            // multiply the documented constant by 10^5 because we want the pressures
            // in [Pa], not in [bar]
            const Scalar Uconst = 0.318316 * 1e5;

            // compute the oil-water Leverett factor.
            const auto& jfuncFlag = jfunc.flag();
            if (jfuncFlag == Opm::JFunc::Flag::WATER || jfuncFlag == Opm::JFunc::Flag::BOTH) {
                // note that we use the surface tension in terms of [dyn/cm]
                Scalar gamma =
                    jfunc.owSurfaceTension();
                pcowLeverettFactor = commonFactor*gamma*Uconst;
            }

            // compute the gas-oil Leverett factor.
            if (jfuncFlag == Opm::JFunc::Flag::GAS || jfuncFlag == Opm::JFunc::Flag::BOTH) {
                // note that we use the surface tension in terms of [dyn/cm]
                Scalar gamma =
                    jfunc.goSurfaceTension();
                pcgoLeverettFactor = commonFactor*gamma*Uconst;
            }
        }
    }
#endif

private:
#if HAVE_ECL_INPUT
    void extractUnscaledSgof_(const Opm::SgofTable& sgofTable)
    {
        // minimum gas and oil-in-gas-oil saturation
        Sgl = sgofTable.getSgColumn().front();
        Sogl = 1.0 - sgofTable.getSgColumn().back();

        // maximum gas and oil-in-gas-oil saturation
        Sgu = sgofTable.getSgColumn().back();
        Sogu = 1.0 - sgofTable.getSgColumn().front();

        // critical gas saturation
        Sgcr = 0.0;
        for (size_t rowIdx = 0; rowIdx < sgofTable.numRows(); ++ rowIdx) {
            if (sgofTable.getKrgColumn()[rowIdx] > krCriticalEps)
                break;

            Sgcr = sgofTable.getSgColumn()[rowIdx];
        }

        // critical oil saturation of gas-oil system
        Sogcr = 0.0;
        for (int rowIdx = static_cast<int>(sgofTable.numRows() - 1);
             rowIdx >= 0;
             -- rowIdx)
        {
            if (sgofTable.getKrogColumn()[static_cast<size_t>(rowIdx)] > krCriticalEps)
                break;

            Sogcr = 1.0 - sgofTable.getSgColumn()[static_cast<size_t>(rowIdx)];
        }

        // maximum gas-oil capillary pressure
        maxPcgo = sgofTable.getPcogColumn().back();

        // maximum gas-* relperms
        maxKrg = sgofTable.getKrgColumn().back();
        maxKrog = sgofTable.getKrogColumn().front();
    }

    void extractUnscaledSlgof_(const Opm::SlgofTable& slgofTable)
    {
        // minimum gas and oil-in-gas-oil saturation
        Sgl = 1.0 - slgofTable.getSlColumn().back();
        Sogl = slgofTable.getSlColumn().front();

        // maximum gas and oil-in-gas-oil saturation
        Sgu = 1.0 - slgofTable.getSlColumn().front();
        Sogu = slgofTable.getSlColumn().back();

        // critical gas saturation
        Sgcr = 0.0;
        for (int rowIdx = static_cast<int>(slgofTable.numRows()) - 1;
             rowIdx >= 0;
             -- rowIdx)
        {
            if (slgofTable.getKrgColumn()[static_cast<size_t>(rowIdx)] > krCriticalEps)
                break;

            Sgcr = 1 - slgofTable.getSlColumn()[static_cast<size_t>(rowIdx)];
        }

        // critical oil saturation of gas-oil system
        Sogcr = 0.0;
        for (size_t rowIdx = 0; rowIdx < slgofTable.numRows(); ++ rowIdx) {
            if (slgofTable.getKrogColumn()[rowIdx] > krCriticalEps)
                break;

            Sogcr = slgofTable.getSlColumn()[rowIdx];
        }

        // maximum gas-oil capillary pressure
        maxPcgo = slgofTable.getPcogColumn().front();

        // maximum gas-* relperms
        maxKrg = slgofTable.getKrgColumn().front();
        maxKrog = slgofTable.getKrogColumn().back();
    }

    void extractUnscaledSwof_(const Opm::SwofTable& swofTable)
    {
        // connate saturations
        Swl = swofTable.getSwColumn().front();
        Sowl = 1.0 - swofTable.getSwColumn().back();

        // maximum water and oil-in-oil-water saturations
        Swu = swofTable.getSwColumn().back();
        Sowu = 1.0 - swofTable.getSwColumn().front();

        // critical water saturation
        Swcr = 0.0;
        for (size_t rowIdx = 0; rowIdx < swofTable.numRows(); ++ rowIdx) {
            if (swofTable.getKrwColumn()[rowIdx] > krCriticalEps)
                break;

            Swcr = swofTable.getSwColumn()[rowIdx];
        }

        // critical oil saturation of oil-water system
        Sowcr = 0.0;
        for (int rowIdx = static_cast<int>(swofTable.numRows()) - 1;
             rowIdx >= 0;
             -- rowIdx)
        {
            if (swofTable.getKrowColumn()[static_cast<size_t>(rowIdx)] > krCriticalEps)
                break;

            Sowcr = 1.0 - swofTable.getSwColumn()[static_cast<size_t>(rowIdx)];
        }

        // maximum oil-water capillary pressures
        maxPcow = swofTable.getPcowColumn().front();

        // maximum water-* relative permeabilities
        maxKrw = swofTable.getKrwColumn().back();
        maxKrow = swofTable.getKrowColumn().front();
    }

    void extractUnscaledSwfn_(const Opm::SwfnTable& swfnTable)
    {
        // connate water saturation
        Swl = swfnTable.getSwColumn().front();

        // maximum water saturation
        Swu = swfnTable.getSwColumn().back();

        // critical water saturation
        Swcr = 0.0;
        for (size_t rowIdx = 0; rowIdx < swfnTable.numRows(); ++ rowIdx) {
            if (swfnTable.getKrwColumn()[rowIdx] > krCriticalEps)
                break;

            Swcr = swfnTable.getSwColumn()[rowIdx];
        }

        // maximum oil-water capillary pressure
        maxPcow = swfnTable.getPcowColumn().front();

        // maximum water relative permeability
        maxKrw = swfnTable.getKrwColumn().back();
    }

    void extractUnscaledSgfn_(const Opm::SgfnTable& sgfnTable)
    {
        // connate gas saturation
        Sgl = sgfnTable.getSgColumn().front();

        // maximum gas saturations
        Sgu = sgfnTable.getSgColumn().back();
        Sogu = 1 - sgfnTable.getSgColumn().front();

        // critical gas saturation
        Sgcr = 0.0;
        for (size_t rowIdx = 0; rowIdx < sgfnTable.numRows(); ++ rowIdx) {
            if (sgfnTable.getKrgColumn()[rowIdx] > krCriticalEps)
                break;

            Sgcr = sgfnTable.getSgColumn()[rowIdx];
        }

        // maximum capillary pressure
        maxPcgo = sgfnTable.getPcogColumn().back();

        // maximum relative gas permeability
        maxKrg = sgfnTable.getKrgColumn().back();
    }

    void extractUnscaledSof3_(const Opm::Sof3Table& sof3Table)
    {
        // connate oil saturations
        Sowl = sof3Table.getSoColumn().front() + Sgl;
        Sogl = sof3Table.getSoColumn().front() + Swl;

        // maximum oil saturations
        Sowu = sof3Table.getSoColumn().back();

        // critical oil saturation of oil-water system
        Sowcr = 0.0;
        for (size_t rowIdx = 0 ; rowIdx < sof3Table.numRows(); ++ rowIdx) {
            if (sof3Table.getKrowColumn()[rowIdx] > krCriticalEps) {
                break;
            }

            Sowcr = sof3Table.getSoColumn()[rowIdx];
        }

        // critical oil saturation of gas-oil system
        Sogcr = 0.0;
        for (size_t rowIdx = 0 ; rowIdx < sof3Table.numRows(); ++ rowIdx) {
            if (sof3Table.getKrogColumn()[rowIdx] > krCriticalEps)
                break;

            Sogcr = sof3Table.getSoColumn()[rowIdx];
        }

        // maximum relative oil permeabilities
        maxKrow = sof3Table.getKrowColumn().back();
        maxKrog = sof3Table.getKrogColumn().back();
    }

    void extractUnscaledSof2_(const Opm::Sof2Table& sof2Table)
    {
        // connate oil saturations
        Sowl = sof2Table.getSoColumn().front() + Sgl;
        Sogl = sof2Table.getSoColumn().front() + Swl;

        // maximum oil saturations
        Sowu = sof2Table.getSoColumn().back();

        // critical oil saturation of oil-water system or critical oil saturation of
        // gas-oil system
        Sowcr = 0.0;
        for (size_t rowIdx = 0 ; rowIdx < sof2Table.numRows(); ++ rowIdx) {
            if (sof2Table.getKroColumn()[rowIdx] > krCriticalEps) {
                break;
            }

            Sowcr = sof2Table.getSoColumn()[rowIdx];
        }
        Sogcr = Sowcr;

        // maximum relative oil permeabilities
        maxKrow = sof2Table.getKroColumn().back();
        maxKrog = maxKrow;
    }
#endif // HAVE_ECL_INPUT

    void extractGridPropertyValue_(Scalar& targetValue,
                                   const std::vector<double>* propData,
                                   unsigned cartesianCellIdx)
    {
        if (!propData)
            return;

        targetValue = (*propData)[cartesianCellIdx];
    }
};

/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Represents the points on the X and Y axis to be scaled if endpoint scaling is
 *        used.
 */
template <class Scalar>
class EclEpsScalingPoints
{
public:
    /*!
     * \brief Assigns the scaling points which actually ought to be used.
     */
    void init(const EclEpsScalingPointsInfo<Scalar>& epsInfo,
              const EclEpsConfig& config,
              EclTwoPhaseSystemType epsSystemType)
    {
        if (epsSystemType == EclOilWaterSystem) {
            // saturation scaling for capillary pressure
            saturationPcPoints_[0] = epsInfo.Swl;
            saturationPcPoints_[1] = epsInfo.Swu;

            // krw saturation scaling endpoints
            if (config.enableThreePointKrSatScaling()) {
                saturationKrwPoints_[0] = epsInfo.Swcr;
                saturationKrwPoints_[1] = 1.0 - epsInfo.Sowcr - epsInfo.Sgl;
                saturationKrwPoints_[2] = epsInfo.Swu;
            }
            else {
                saturationKrwPoints_[0] = epsInfo.Swcr;
                saturationKrwPoints_[1] = epsInfo.Swu;
            }

            // krn saturation scaling endpoints (with the non-wetting phase being oil).
            // because opm-material specifies non-wetting phase relperms in terms of the
            // wetting phase saturations, the code here uses 1 minus the values specified
            // by the Eclipse TD and the order of the scaling points is reversed
            if (config.enableThreePointKrSatScaling()) {
                saturationKrnPoints_[2] = 1.0 - epsInfo.Sowcr;
                saturationKrnPoints_[1] = epsInfo.Swcr + epsInfo.Sgl;
                saturationKrnPoints_[0] = epsInfo.Swl + epsInfo.Sgl;
            }
            else {
                saturationKrnPoints_[1] = 1 - epsInfo.Sowcr;
                saturationKrnPoints_[0] = epsInfo.Swl + epsInfo.Sgl;
            }

            if (config.enableLeverettScaling())
                maxPcnwOrLeverettFactor_ = epsInfo.pcowLeverettFactor;
            else
                maxPcnwOrLeverettFactor_ = epsInfo.maxPcow;
            maxKrw_ = epsInfo.maxKrw;
            maxKrn_ = epsInfo.maxKrow;
        }
        else {
            assert(epsSystemType == EclGasOilSystem);

            // saturation scaling for capillary pressure
            saturationPcPoints_[0] = 1.0 - epsInfo.Sgu;
            saturationPcPoints_[1] = 1.0 - epsInfo.Sgl;

            // krw saturation scaling endpoints
            if (config.enableThreePointKrSatScaling()) {
                saturationKrwPoints_[0] = epsInfo.Sogcr;
                saturationKrwPoints_[1] = 1 - epsInfo.Sgcr - epsInfo.Swl;
                saturationKrwPoints_[2] = 1 - epsInfo.Swl - epsInfo.Sgl;
            }
            else {
                saturationKrwPoints_[0] = epsInfo.Sogcr;
                saturationKrwPoints_[1] = 1 - epsInfo.Swl - epsInfo.Sgl;
            }

            // krn saturation scaling endpoints (with the non-wetting phase being gas).
            // because opm-material specifies non-wetting phase relperms in terms of the
            // wetting phase saturations, the code here uses 1 minus the values specified
            // by the Eclipse TD and the order of the scaling points is reversed
            if (config.enableThreePointKrSatScaling()) {
                saturationKrnPoints_[2] = 1.0 - epsInfo.Sgcr;
                saturationKrnPoints_[1] = epsInfo.Sogcr + epsInfo.Swl;
                saturationKrnPoints_[0] = 1.0 - epsInfo.Sgu;
            }
            else {
                saturationKrnPoints_[1] = 1.0 - epsInfo.Sgcr;
                saturationKrnPoints_[0] = 1.0 - epsInfo.Sgu;
            }

            if (config.enableLeverettScaling())
                maxPcnwOrLeverettFactor_ = epsInfo.pcgoLeverettFactor;
            else
                maxPcnwOrLeverettFactor_ = epsInfo.maxPcgo;

            maxKrw_ = epsInfo.maxKrog;
            maxKrn_ = epsInfo.maxKrg;
        }
    }

    /*!
     * \brief Sets an saturation value for capillary pressure saturation scaling
     */
    void setSaturationPcPoint(unsigned pointIdx, Scalar value)
    { saturationPcPoints_[pointIdx] = value; }

    /*!
     * \brief Returns the points used for capillary pressure saturation scaling
     */
    const std::array<Scalar, 2>& saturationPcPoints() const
    { return saturationPcPoints_; }

    /*!
     * \brief Sets an saturation value for wetting-phase relperm saturation scaling
     */
    void setSaturationKrwPoint(unsigned pointIdx, Scalar value)
    { saturationKrwPoints_[pointIdx] = value; }

    /*!
     * \brief Returns the points used for wetting phase relperm saturation scaling
     */
    const std::array<Scalar, 3>& saturationKrwPoints() const
    { return saturationKrwPoints_; }

    /*!
     * \brief Sets an saturation value for non-wetting phase relperm saturation scaling
     */
    void setSaturationKrnPoint(unsigned pointIdx, Scalar value)
    { saturationKrnPoints_[pointIdx] = value; }

    /*!
     * \brief Returns the points used for non-wetting phase relperm saturation scaling
     */
    const std::array<Scalar, 3>& saturationKrnPoints() const
    { return saturationKrnPoints_; }

    /*!
     * \brief Sets the maximum capillary pressure
     */
    void setMaxPcnw(Scalar value)
    { maxPcnwOrLeverettFactor_ = value; }

    /*!
     * \brief Returns the maximum capillary pressure
     */
    Scalar maxPcnw() const
    { return maxPcnwOrLeverettFactor_; }

    /*!
     * \brief Sets the Leverett scaling factor for capillary pressure
     */
    void setLeverettFactor(Scalar value)
    { maxPcnwOrLeverettFactor_ = value; }

    /*!
     * \brief Returns the Leverett scaling factor for capillary pressure
     */
    Scalar leverettFactor() const
    { return maxPcnwOrLeverettFactor_; }

    /*!
     * \brief Sets the maximum wetting phase relative permeability
     */
    void setMaxKrw(Scalar value)
    { maxKrw_ = value; }

    /*!
     * \brief Returns the maximum wetting phase relative permeability
     */
    Scalar maxKrw() const
    { return maxKrw_; }

    /*!
     * \brief Sets the maximum wetting phase relative permeability
     */
    void setMaxKrn(Scalar value)
    { maxKrn_ = value; }

    /*!
     * \brief Returns the maximum wetting phase relative permeability
     */
    Scalar maxKrn() const
    { return maxKrn_; }

    void print() const
    {
        std::cout << "    saturationKrnPoints_[0]: " << saturationKrnPoints_[0] << "\n"
                  << "    saturationKrnPoints_[1]: " << saturationKrnPoints_[1] << "\n"
                  << "    saturationKrnPoints_[2]: " << saturationKrnPoints_[2] << "\n";
    }

private:
    // The the points used for the "y-axis" scaling of capillary pressure
    Scalar maxPcnwOrLeverettFactor_;

    // The the points used for the "y-axis" scaling of wetting phase relative permability
    Scalar maxKrw_;

    // The the points used for the "y-axis" scaling of non-wetting phase relative permability
    Scalar maxKrn_;

    // The the points used for saturation ("x-axis") scaling of capillary pressure
    std::array<Scalar, 2> saturationPcPoints_;

    // The the points used for saturation ("x-axis") scaling of wetting phase relative permeability
    std::array<Scalar, 3> saturationKrwPoints_;

    // The the points used for saturation ("x-axis") scaling of non-wetting phase relative permeability
    std::array<Scalar, 3> saturationKrnPoints_;
};

} // namespace Opm

#endif
