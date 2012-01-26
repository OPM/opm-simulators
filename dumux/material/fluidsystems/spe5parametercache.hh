/*****************************************************************************
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 *
 * \brief Specifies the parameters required by the SPE5 problem which
 *        are dependend on the thermodynamic state.
 */
#ifndef SPE5_PARAMETER_CACHE_HH
#define SPE5_PARAMETER_CACHE_HH

#include <dumux/material/components/h2o.hh>
#include <dumux/material/fluidsystems/parametercachebase.hh>

#include <dumux/material/eos/pengrobinson.hh>
#include <dumux/material/eos/pengrobinsonparamsmixture.hh>

#include <assert.h>

namespace Dumux
{
/*!
 * \brief Specifies the parameters required by the SPE5 problem which
 *        are dependend on the thermodynamic state.
 */
template <class Scalar, class FluidSystem>
class Spe5ParameterCache
    : public Dumux::ParameterCacheBase<Spe5ParameterCache<Scalar, FluidSystem> > 
{
    typedef Spe5ParameterCache<Scalar, FluidSystem> ThisType;
    typedef Dumux::ParameterCacheBase<ThisType> ParentType;
  
    typedef Dumux::PengRobinson<Scalar> PengRobinson;

    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    enum { wPhaseIdx = FluidSystem::wPhaseIdx };
    enum { oPhaseIdx = FluidSystem::oPhaseIdx };
    enum { gPhaseIdx = FluidSystem::gPhaseIdx };

public:
    // types of the parameter objects for each phase
    typedef Dumux::PengRobinsonParamsMixture<Scalar, FluidSystem, oPhaseIdx, /*useSpe5=*/true> OilPhaseParams;
    typedef Dumux::PengRobinsonParamsMixture<Scalar, FluidSystem, gPhaseIdx, /*useSpe5=*/true> GasPhaseParams;

    /*!
     * \brief The constructor
     */
    Spe5ParameterCache()
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            VmUpToDate_[phaseIdx] = false;
            Valgrind::SetUndefined(Vm_[phaseIdx]);
        }
    };

    /*!
     * \brief Update all parameters required by the fluid system to
     *        calculate some quantities for the phase.
     */
    template <class FluidState>
    void updatePhase(const FluidState &fs, 
                     int phaseIdx, 
                     int except = ParentType::None)
    { 
        updateEosParams(fs, phaseIdx, except);

        // if we don't need to recalculate the molar volume, we exit
        // here
        if (VmUpToDate_[phaseIdx])
            return;

        // update the phase's molar volume
        updateMolarVolume_(fs, phaseIdx);
    }
    
    /*!
     * \brief Update all cached parameters of a specific fluid phase
     *        which depend on the mole fraction of a single component
     *
     * *Only* use this method if just a single component's
     * concentration changed between two update*() calls. If more than
     * one concentration changed, call updatePhaseComposition() of
     * updatePhase()!
     */
    template <class FluidState>
    void updateSingleMoleFraction(const FluidState &fs,
                                  int phaseIdx,
                                  int compIdx)
    {
        if (phaseIdx == oPhaseIdx)
            oilPhaseParams_.updateSingleMoleFraction(fs,
                                                     compIdx,
                                                     fs.moleFraction(phaseIdx, compIdx)
                                                     - moleFrac_[phaseIdx][compIdx]);
        else if (phaseIdx == gPhaseIdx)
            gasPhaseParams_.updateSingleMoleFraction(fs,
                                                     compIdx,
                                                     fs.moleFraction(phaseIdx, compIdx)
                                                     - moleFrac_[phaseIdx][compIdx]);

        // update the mole fraction which the parameters are
        // calculated for
        moleFrac_[phaseIdx][compIdx] = fs.moleFraction(phaseIdx, compIdx);
        
        // update the phase's molar volume
        updateMolarVolume_(fs, phaseIdx);
    }

    /*!
     * \brief The Peng-Robinson attractive parameter for a phase.
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
     */
    template <class FluidState>
    void updateEosParams(const FluidState &fs, 
                         int phaseIdx,
                         int exceptQuantities = ParentType::None)
    { 
        if (!(exceptQuantities & ParentType::Temperature))
        {
            updatePure_(fs, phaseIdx);
            updateMix_(fs, phaseIdx);

            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                moleFrac_[phaseIdx][compIdx] = fs.moleFraction(phaseIdx, compIdx);
            VmUpToDate_[phaseIdx] = false;
        }
        else if (!(exceptQuantities & ParentType::Composition))
        {
            updateMix_(fs, phaseIdx);
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                moleFrac_[phaseIdx][compIdx] = fs.moleFraction(phaseIdx, compIdx);
            VmUpToDate_[phaseIdx] = false;
        }
        else if (!(exceptQuantities & ParentType::Pressure))
            VmUpToDate_[phaseIdx] = false;
    }

protected:
    /*!
     * \brief Update all parameters of a phase which only depend on 
     *        temperature and/or pressure.
     *
     * This usually means the parameters for the pure components.
     */
    template <class FluidState>
    void updatePure_(const FluidState &fs, int phaseIdx)
    {
        Scalar T = fs.temperature(phaseIdx);
        Scalar p = fs.pressure(phaseIdx);

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
    void updateMix_(const FluidState &fs, int phaseIdx)
    {
        Valgrind::CheckDefined(fs.averageMolarMass(phaseIdx));
        switch (phaseIdx)
        {
        case oPhaseIdx: 
            oilPhaseParams_.updateMix(fs);
            break;
        case gPhaseIdx:
            gasPhaseParams_.updateMix(fs); 
            break;
        case wPhaseIdx:
            break;
        }
    }

    template <class FluidState>
    void updateMolarVolume_(const FluidState &fs, 
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
                PengRobinson::computeMolarVolume(fs,
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
                PengRobinson::computeMolarVolume(fs,
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
            Scalar overPressure = fs.pressure(wPhaseIdx) - 1.013e5; // [Pa]
            Scalar waterDensity = 
                stockTankWaterDensity * (1 + 3.3e-6*overPressure/6894.7573);

            // convert water density [kg/m^3] to molar volume [m^3/mol]
            Vm_[wPhaseIdx] = fs.averageMolarMass(wPhaseIdx)/waterDensity;
        };
        };
    }

    bool VmUpToDate_[numPhases];
    Scalar Vm_[numPhases];
    Scalar moleFrac_[numPhases][numComponents];

    OilPhaseParams oilPhaseParams_;
    GasPhaseParams gasPhaseParams_;
};

} // end namespace Dumux

#endif
