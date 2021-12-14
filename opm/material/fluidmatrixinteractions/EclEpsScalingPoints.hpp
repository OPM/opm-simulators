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
#include "EclEpsGridProperties.hpp"

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Deck/DeckRecord.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>
#include <opm/input/eclipse/EclipseState/Grid/SatfuncPropertyInitializers.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#endif

#include <opm/material/common/Means.hpp>

#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include <algorithm>

namespace Opm {

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
    Scalar Swl; // water
    Scalar Sgl; // gas

    // critical saturations
    Scalar Swcr; // water
    Scalar Sgcr; // gas
    Scalar Sowcr; // oil for the oil-water system
    Scalar Sogcr; // oil for the gas-oil system

    // maximum saturations
    Scalar Swu; // water
    Scalar Sgu; // gas

    // maximum capillary pressures
    Scalar maxPcow; // maximum capillary pressure of the oil-water system
    Scalar maxPcgo; // maximum capillary pressure of the gas-oil system

    // the Leverett capillary pressure scaling factors. (those only make sense for the
    // scaled points, for the unscaled ones they are 1.0.)
    Scalar pcowLeverettFactor;
    Scalar pcgoLeverettFactor;

    // Scaled relative permeabilities at residual displacing saturation
    Scalar Krwr;  // water
    Scalar Krgr;  // gas
    Scalar Krorw; // oil in water-oil system
    Scalar Krorg; // oil in gas-oil system

    // maximum relative permabilities
    Scalar maxKrw; // maximum relative permability of water
    Scalar maxKrow; // maximum relative permability of oil in the oil-water system
    Scalar maxKrog; // maximum relative permability of oil in the gas-oil system
    Scalar maxKrg; // maximum relative permability of gas

    bool operator==(const EclEpsScalingPointsInfo<Scalar>& data) const
    {
        return Swl == data.Swl &&
               Sgl == data.Sgl &&
               Swcr == data.Swcr &&
               Sgcr == data.Sgcr &&
               Sowcr == data.Sowcr &&
               Sogcr == data.Sogcr &&
               Swu == data.Swu &&
               Sgu == data.Sgu &&
               maxPcow == data.maxPcow &&
               maxPcgo == data.maxPcgo &&
               pcowLeverettFactor == data.pcowLeverettFactor &&
               pcgoLeverettFactor == data.pcgoLeverettFactor &&
               Krwr == data.Krwr &&
               Krgr == data.Krgr &&
               Krorw == data.Krorw &&
               Krorg == data.Krorg &&
               maxKrw == data.maxKrw &&
               maxKrow == data.maxKrow &&
               maxKrog == data.maxKrog &&
               maxKrg == data.maxKrg;
    }

    void print() const
    {
        std::cout << "    Swl: " << Swl << '\n'
                  << "    Sgl: " << Sgl << '\n'
                  << "    Swcr: " << Swcr << '\n'
                  << "    Sgcr: " << Sgcr << '\n'
                  << "    Sowcr: " << Sowcr << '\n'
                  << "    Sogcr: " << Sogcr << '\n'
                  << "    Swu: " << Swu << '\n'
                  << "    Sgu: " << Sgu << '\n'
                  << "    maxPcow: " << maxPcow << '\n'
                  << "    maxPcgo: " << maxPcgo << '\n'
                  << "    pcowLeverettFactor: " << pcowLeverettFactor << '\n'
                  << "    pcgoLeverettFactor: " << pcgoLeverettFactor << '\n'
                  << "    Krwr: " << Krwr << '\n'
                  << "    Krgr: " << Krgr << '\n'
                  << "    Krorw: " << Krorw << '\n'
                  << "    Krorg: " << Krorg << '\n'
                  << "    maxKrw: " << maxKrw << '\n'
                  << "    maxKrg: " << maxKrg << '\n'
                  << "    maxKrow: " << maxKrow << '\n'
                  << "    maxKrog: " << maxKrog << '\n';
    }

#if HAVE_ECL_INPUT
    /*!
     * \brief Extract the values of the unscaled scaling parameters.
     *
     * I.e., the values which are used for the nested Fluid-Matrix interactions and which
     * are produced by them.
     */
    void extractUnscaled(const satfunc::RawTableEndPoints& rtep,
                         const satfunc::RawFunctionValues& rfunc,
                         const std::vector<double>::size_type   satRegionIdx)
    {
        this->Swl = rtep.connate.water[satRegionIdx];
        this->Sgl = rtep.connate.gas  [satRegionIdx];

        this->Swcr  = rtep.critical.water       [satRegionIdx];
        this->Sgcr  = rtep.critical.gas         [satRegionIdx];
        this->Sowcr = rtep.critical.oil_in_water[satRegionIdx];
        this->Sogcr = rtep.critical.oil_in_gas  [satRegionIdx];

        this->Swu = rtep.maximum.water[satRegionIdx];
        this->Sgu = rtep.maximum.gas  [satRegionIdx];

        this->maxPcgo = rfunc.pc.g[satRegionIdx];
        this->maxPcow = rfunc.pc.w[satRegionIdx];

        // there are no "unscaled" Leverett factors, so we just set them to 1.0
        this->pcowLeverettFactor = 1.0;
        this->pcgoLeverettFactor = 1.0;

        this->Krwr    = rfunc.krw.r [satRegionIdx];
        this->Krgr    = rfunc.krg.r [satRegionIdx];
        this->Krorw   = rfunc.kro.rw[satRegionIdx];
        this->Krorg   = rfunc.kro.rg[satRegionIdx];

        this->maxKrw  = rfunc.krw.max[satRegionIdx];
        this->maxKrow = rfunc.kro.max[satRegionIdx];
        this->maxKrog = rfunc.kro.max[satRegionIdx];
        this->maxKrg  = rfunc.krg.max[satRegionIdx];
    }

    void update(Scalar& targetValue, const double * value_ptr) {
        if (value_ptr)
            targetValue = *value_ptr;
    }

    /*!
     * \brief Extract the values of the scaled scaling parameters.
     *
     * I.e., the values which are "seen" by the physical model.
     */
    void extractScaled(const EclipseState& eclState,
                       const EclEpsGridProperties& epsProperties,
                       unsigned activeIndex)
    {
        // overwrite the unscaled values with the values for the cell if it is
        // explicitly specified by the corresponding keyword.
        update(Swl,     epsProperties.swl(activeIndex));
        update(Sgl,     epsProperties.sgl(activeIndex));
        update(Swcr,    epsProperties.swcr(activeIndex));
        update(Sgcr,    epsProperties.sgcr(activeIndex));

        update(Sowcr,   epsProperties.sowcr(activeIndex));
        update(Sogcr,   epsProperties.sogcr(activeIndex));
        update(Swu,     epsProperties.swu(activeIndex));
        update(Sgu,     epsProperties.sgu(activeIndex));
        update(maxPcow, epsProperties.pcw(activeIndex));
        update(maxPcgo, epsProperties.pcg(activeIndex));

        update(this->Krwr,  epsProperties.krwr(activeIndex));
        update(this->Krgr,  epsProperties.krgr(activeIndex));
        update(this->Krorw, epsProperties.krorw(activeIndex));
        update(this->Krorg, epsProperties.krorg(activeIndex));

        update(maxKrw,  epsProperties.krw(activeIndex));
        update(maxKrg,  epsProperties.krg(activeIndex));
        update(maxKrow, epsProperties.kro(activeIndex));
        update(maxKrog, epsProperties.kro(activeIndex));

        // compute the Leverett capillary pressure scaling factors if applicable.  note
        // that this needs to be done using non-SI units to make it correspond to the
        // documentation.
        pcowLeverettFactor = 1.0;
        pcgoLeverettFactor = 1.0;
        if (eclState.getTableManager().useJFunc()) {
            const auto& jfunc = eclState.getTableManager().getJFunc();
            const auto& jfuncDir = jfunc.direction();

            Scalar perm;
            if (jfuncDir == JFunc::Direction::X)
                perm = epsProperties.permx(activeIndex);
            else if (jfuncDir == JFunc::Direction::Y)
                perm = epsProperties.permy(activeIndex);
            else if (jfuncDir == JFunc::Direction::Z)
                perm = epsProperties.permz(activeIndex);
            else if (jfuncDir == JFunc::Direction::XY)
                // TODO: verify that this really is the arithmetic mean. (the
                // documentation just says that the "average" should be used, IMO the
                // harmonic mean would be more appropriate because that's what's usually
                // applied when calculating the fluxes.)
            {
                double permx = epsProperties.permx(activeIndex);
                double permy = epsProperties.permy(activeIndex);
                perm = arithmeticMean(permx, permy);
            } else
                throw std::runtime_error("Illegal direction indicator for the JFUNC "
                                         "keyword ("+std::to_string(int(jfuncDir))+")");

            // convert permeability from m^2 to mD
            perm *= 1.01325e15;

            Scalar poro = epsProperties.poro(activeIndex);
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
            if (jfuncFlag == JFunc::Flag::WATER || jfuncFlag == JFunc::Flag::BOTH) {
                // note that we use the surface tension in terms of [dyn/cm]
                Scalar gamma =
                    jfunc.owSurfaceTension();
                pcowLeverettFactor = commonFactor*gamma*Uconst;
            }

            // compute the gas-oil Leverett factor.
            if (jfuncFlag == JFunc::Flag::GAS || jfuncFlag == JFunc::Flag::BOTH) {
                // note that we use the surface tension in terms of [dyn/cm]
                Scalar gamma =
                    jfunc.goSurfaceTension();
                pcgoLeverettFactor = commonFactor*gamma*Uconst;
            }
        }
    }
#endif

private:
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
            saturationPcPoints_[2] = saturationPcPoints_[1] = epsInfo.Swu;

            // krw saturation scaling endpoints
            saturationKrwPoints_[0] = epsInfo.Swcr;
            saturationKrwPoints_[1] = 1.0 - epsInfo.Sowcr - epsInfo.Sgl;
            saturationKrwPoints_[2] = epsInfo.Swu;

            // krn saturation scaling endpoints (with the non-wetting phase being oil).
            // because opm-material specifies non-wetting phase relperms in terms of the
            // wetting phase saturations, the code here uses 1 minus the values specified
            // by the Eclipse TD and the order of the scaling points is reversed
            saturationKrnPoints_[2] = 1.0 - epsInfo.Sowcr;
            saturationKrnPoints_[1] = epsInfo.Swcr + epsInfo.Sgl;
            saturationKrnPoints_[0] = epsInfo.Swl + epsInfo.Sgl;

            if (config.enableLeverettScaling())
                maxPcnwOrLeverettFactor_ = epsInfo.pcowLeverettFactor;
            else
                maxPcnwOrLeverettFactor_ = epsInfo.maxPcow;

            Krwr_   = epsInfo.Krwr;
            Krnr_   = epsInfo.Krorw;

            maxKrw_ = epsInfo.maxKrw;
            maxKrn_ = epsInfo.maxKrow;
        }
        else {
            assert((epsSystemType == EclGasOilSystem) ||(epsSystemType == EclGasWaterSystem) );

            // saturation scaling for capillary pressure
            saturationPcPoints_[0] = 1.0 - epsInfo.Swl - epsInfo.Sgu;
            saturationPcPoints_[2] = saturationPcPoints_[1] = 1.0 - epsInfo.Swl - epsInfo.Sgl;

            // krw saturation scaling endpoints
            saturationKrwPoints_[0] = epsInfo.Sogcr;
            saturationKrwPoints_[1] = 1.0 - epsInfo.Sgcr - epsInfo.Swl;
            saturationKrwPoints_[2] = 1.0 - epsInfo.Swl - epsInfo.Sgl;

            // krn saturation scaling endpoints (with the non-wetting phase being gas).
            //
            // As opm-material specifies non-wetting phase relative
            // permeabilities in terms of the wetting phase saturations, the
            // code here uses (1-SWL) minus the values specified by the
            // ECLIPSE TD and the order of the scaling points is reversed.
            saturationKrnPoints_[2] = 1.0 - epsInfo.Swl - epsInfo.Sgcr;
            saturationKrnPoints_[1] = epsInfo.Sogcr;
            saturationKrnPoints_[0] = 1.0 - epsInfo.Swl - epsInfo.Sgu;

            if (config.enableLeverettScaling())
                maxPcnwOrLeverettFactor_ = epsInfo.pcgoLeverettFactor;
            else
                maxPcnwOrLeverettFactor_ = epsInfo.maxPcgo;

            Krwr_   = epsInfo.Krorg;
            Krnr_   = epsInfo.Krgr;

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
    const std::array<Scalar, 3>& saturationPcPoints() const
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
     * \brief Set wetting-phase relative permeability at residual saturation
     * of non-wetting phase.
     */
    void setKrwr(Scalar value)
    { this->Krwr_ = value; }

    /*!
     * \brief Returns wetting-phase relative permeability at residual
     * saturation of non-wetting phase.
     */
    Scalar krwr() const
    { return this->Krwr_; }

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
     * \brief Set non-wetting phase relative permeability at residual
     * saturation of wetting phase.
     */
    void setKrnr(Scalar value)
    { this->Krnr_ = value; }

    /*!
     * \brief Returns non-wetting phase relative permeability at residual
     * saturation of wetting phase.
     */
    Scalar krnr() const
    { return this->Krnr_; }

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
    // Points used for vertical scaling of capillary pressure
    Scalar maxPcnwOrLeverettFactor_;

    // Maximum wetting phase relative permability value.
    Scalar maxKrw_;

    // Scaled wetting phase relative permeability value at residual
    // saturation of non-wetting phase.
    Scalar Krwr_;

    // Maximum non-wetting phase relative permability value
    Scalar maxKrn_;

    // Scaled non-wetting phase relative permeability value at residual
    // saturation of wetting phase.
    Scalar Krnr_;

    // The the points used for saturation ("x-axis") scaling of capillary pressure
    std::array<Scalar, 3> saturationPcPoints_;

    // The the points used for saturation ("x-axis") scaling of wetting phase relative permeability
    std::array<Scalar, 3> saturationKrwPoints_;

    // The the points used for saturation ("x-axis") scaling of non-wetting phase relative permeability
    std::array<Scalar, 3> saturationKrnPoints_;
};

} // namespace Opm

#endif
