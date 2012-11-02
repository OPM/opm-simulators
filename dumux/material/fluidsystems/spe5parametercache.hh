// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \copydoc Dumux::Spe5ParameterCache
 */
#ifndef DUMUX_SPE5_PARAMETER_CACHE_HH
#define DUMUX_SPE5_PARAMETER_CACHE_HH

#include <cassert>

#include <dumux/material/components/h2o.hh>
#include <dumux/material/fluidsystems/parametercachebase.hh>

#include <dumux/material/eos/pengrobinson.hh>
#include <dumux/material/eos/pengrobinsonparamsmixture.hh>

namespace Dumux {

/*!
 * \ingroup Fluidsystems
 * \brief Specifies the parameter cache used by the SPE-5 fluid system.
 */
template <class Scalar, class FluidSystem>
class Spe5ParameterCache
    : public Dumux::ParameterCacheBase<Spe5ParameterCache<Scalar, FluidSystem> >
{
    typedef Spe5ParameterCache<Scalar, FluidSystem> ThisType;
    typedef Dumux::ParameterCacheBase<ThisType> ParentType;

    typedef Dumux::PengRobinson<Scalar> PengRobinson;

    enum { numPhases = FluidSystem::numPhases };

    enum { wPhaseIdx = FluidSystem::wPhaseIdx };
    enum { oPhaseIdx = FluidSystem::oPhaseIdx };
    enum { gPhaseIdx = FluidSystem::gPhaseIdx };

public:
    //! The cached parameters for the oil phase
    typedef Dumux::PengRobinsonParamsMixture<Scalar, FluidSystem, oPhaseIdx, /*useSpe5=*/true> OilPhaseParams;
    //! The cached parameters for the gas phase
    typedef Dumux::PengRobinsonParamsMixture<Scalar, FluidSystem, gPhaseIdx, /*useSpe5=*/true> GasPhaseParams;

    Spe5ParameterCache()
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            VmUpToDate_[phaseIdx] = false;
            Valgrind::SetUndefined(Vm_[phaseIdx]);
        }
    }

    //! \copydoc ParameterCacheBase::updatePhase
    template <class FluidState>
    void updatePhase(const FluidState &fluidState,
                     int phaseIdx,
                     int exceptQuantities = ParentType::None)
    {
        updateEosParams(fluidState, phaseIdx, exceptQuantities);

        // if we don't need to recalculate the molar volume, we exit
        // here
        if (VmUpToDate_[phaseIdx])
            return;

        // update the phase's molar volume
        updateMolarVolume_(fluidState, phaseIdx);
    }

    //! \copydoc ParameterCacheBase::updateSingleMoleFraction
    template <class FluidState>
    void updateSingleMoleFraction(const FluidState &fluidState,
                                  int phaseIdx,
                                  int compIdx)
    {
        if (phaseIdx == oPhaseIdx)
            oilPhaseParams_.updateSingleMoleFraction(fluidState, compIdx);
        else if (phaseIdx == gPhaseIdx)
            gasPhaseParams_.updateSingleMoleFraction(fluidState, compIdx);

        // update the phase's molar volume
        updateMolarVolume_(fluidState, phaseIdx);
    }

    /*!
     * \brief The Peng-Robinson attractive parameter for a phase.
     *
     * \param phaseIdx The fluid phase of interest
     */
    Scalar a(int phaseIdx) const
    {
        switch (phaseIdx)
        {
        case oPhaseIdx: return oilPhaseParams_.a();
        case gPhaseIdx: return gasPhaseParams_.a();
        default:
            DUNE_THROW(Dune::InvalidStateException,
                       "The a() parameter is only defined for "
                       "oil and gas phases");
        };
    }

    /*!
     * \brief The Peng-Robinson covolume for a phase.
     *
     * \param phaseIdx The fluid phase of interest
     */
    Scalar b(int phaseIdx) const
    {
        switch (phaseIdx)
        {
        case oPhaseIdx: return oilPhaseParams_.b();
        case gPhaseIdx: return gasPhaseParams_.b();
        default:
            DUNE_THROW(Dune::InvalidStateException,
                       "The b() parameter is only defined for "
                       "oil and gas phases");
        };
    }

    /*!
     * \brief The Peng-Robinson attractive parameter for a pure
     *        component given the same temperature and pressure of the
     *        phase.
     *
     * \param phaseIdx The fluid phase of interest
     * \param compIdx The component phase of interest
     */
    Scalar aPure(int phaseIdx, int compIdx) const
    {
        switch (phaseIdx)
        {
        case oPhaseIdx: return oilPhaseParams_.pureParams(compIdx).a();
        case gPhaseIdx: return gasPhaseParams_.pureParams(compIdx).a();
        default:
            DUNE_THROW(Dune::InvalidStateException,
                       "The a() parameter is only defined for "
                       "oil and gas phases");
        };
    }

    /*!
     * \brief The Peng-Robinson covolume for a pure component given
     *        the same temperature and pressure of the phase.
     *
     * \param phaseIdx The fluid phase of interest
     * \param compIdx The component phase of interest
     */
    Scalar bPure(int phaseIdx, int compIdx) const
    {
        switch (phaseIdx)
        {
        case oPhaseIdx: return oilPhaseParams_.pureParams(compIdx).b();
        case gPhaseIdx: return gasPhaseParams_.pureParams(compIdx).b();
        default:
            DUNE_THROW(Dune::InvalidStateException,
                       "The b() parameter is only defined for "
                       "oil and gas phases");
        };
    }

    /*!
     * \brief Returns the molar volume of a phase [m^3/mol]
     *
     * \param phaseIdx The fluid phase of interest
     */
    Scalar molarVolume(int phaseIdx) const
    { assert(VmUpToDate_[phaseIdx]); return Vm_[phaseIdx]; }


    /*!
     * \brief Returns the Peng-Robinson mixture parameters for the oil
     *        phase.
     */
    const OilPhaseParams &oilPhaseParams() const
    { return oilPhaseParams_; }

    /*!
     * \brief Returns the Peng-Robinson mixture parameters for the gas
     *        phase.
     */
    const GasPhaseParams &gasPhaseParams() const
    { return gasPhaseParams_; }

    /*!
     * \brief Update all parameters required by the equation of state to
     *        calculate some quantities for the phase.
     *
     * \param fluidState The representation of the thermodynamic system of interest.
     * \param phaseIdx The index of the fluid phase of interest.
     * \param exceptQuantities The quantities of the fluid state that have not changed since the last update.
     */
    template <class FluidState>
    void updateEosParams(const FluidState &fluidState,
                         int phaseIdx,
                         int exceptQuantities = ParentType::None)
    {
        if (!(exceptQuantities & ParentType::Temperature))
        {
            updatePure_(fluidState, phaseIdx);
            updateMix_(fluidState, phaseIdx);
            VmUpToDate_[phaseIdx] = false;
        }
        else if (!(exceptQuantities & ParentType::Composition))
        {
            updateMix_(fluidState, phaseIdx);
            VmUpToDate_[phaseIdx] = false;
        }
        else if (!(exceptQuantities & ParentType::Pressure)) {
            VmUpToDate_[phaseIdx] = false;
        }
    }

protected:
    /*!
     * \brief Update all parameters of a phase which only depend on
     *        temperature and/or pressure.
     *
     * This usually means the parameters for the pure components.
     */
    template <class FluidState>
    void updatePure_(const FluidState &fluidState, int phaseIdx)
    {
        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        switch (phaseIdx)
        {
        case oPhaseIdx: oilPhaseParams_.updatePure(T, p); break;
        case gPhaseIdx: gasPhaseParams_.updatePure(T, p); break;
            //case wPhaseIdx: waterPhaseParams_.updatePure(phaseIdx, temperature, pressure);break;
        }
    }

    /*!
     * \brief Update all parameters of a phase which depend on the
     *        fluid composition. It is assumed that updatePure() has
     *        been called before this method.
     *
     * Here, the mixing rule kicks in.
     */
    template <class FluidState>
    void updateMix_(const FluidState &fluidState, int phaseIdx)
    {
        Valgrind::CheckDefined(fluidState.averageMolarMass(phaseIdx));
        switch (phaseIdx)
        {
        case oPhaseIdx:
            oilPhaseParams_.updateMix(fluidState);
            break;
        case gPhaseIdx:
            gasPhaseParams_.updateMix(fluidState);
            break;
        case wPhaseIdx:
            break;
        }
    }

    template <class FluidState>
    void updateMolarVolume_(const FluidState &fluidState,
                            int phaseIdx)
    {
        VmUpToDate_[phaseIdx] = true;

        // calculate molar volume of the phase (we will need this for the
        // fugacity coefficients and the density anyway)
        switch (phaseIdx) {
        case gPhaseIdx: {
            // calculate molar volumes for the given composition. although
            // this isn't a Peng-Robinson parameter strictly speaking, the
            // molar volume appears in basically every quantity the fluid
            // system can get queried, so it is okay to calculate it
            // here...
            Vm_[gPhaseIdx] =
                PengRobinson::computeMolarVolume(fluidState,
                                                 *this,
                                                 phaseIdx,
                                                 /*isGasPhase=*/true);
        }
        case oPhaseIdx: {
            // calculate molar volumes for the given composition. although
            // this isn't a Peng-Robinson parameter strictly speaking, the
            // molar volume appears in basically every quantity the fluid
            // system can get queried, so it is okay to calculate it
            // here...
            Vm_[oPhaseIdx] =
                PengRobinson::computeMolarVolume(fluidState,
                                                 *this,
                                                 phaseIdx,
                                                 /*isGasPhase=*/false);

        }
        case wPhaseIdx: {
            // Density of water in the stock tank (i.e. atmospheric
            // pressure) is specified as 62.4 lb/ft^3 by the SPE-5
            // paper. Also 1 lb = 0.4535923 and 1 ft = 0.3048 m.
            const Scalar stockTankWaterDensity = 62.4 * 0.45359237 / 0.028316847;
            // Water compressibility is specified as 3.3e-6 per psi
            // overpressure, where 1 psi = 6894.7573 Pa
            Scalar overPressure = fluidState.pressure(wPhaseIdx) - 1.013e5; // [Pa]
            Scalar waterDensity =
                stockTankWaterDensity * (1 + 3.3e-6*overPressure/6894.7573);

            // convert water density [kg/m^3] to molar volume [m^3/mol]
            Vm_[wPhaseIdx] = fluidState.averageMolarMass(wPhaseIdx)/waterDensity;
        };
        };
    }

    bool VmUpToDate_[numPhases];
    Scalar Vm_[numPhases];

    OilPhaseParams oilPhaseParams_;
    GasPhaseParams gasPhaseParams_;
};

} // end namespace Dumux

#endif
