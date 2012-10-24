// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 * \brief The mixing rule for the oil and the gas phases of the SPE5 problem.
 *
 * This problem comprises \f$H_2O\f$, \f$C_1\f$, \f$C_3\f$, \f$C_6\f$,
 * \f$C_10\f$, \f$C_15\f$ and \f$C_20\f$ as components.
 *
 * See:
 *
 * R. Reid, et al.: The Properties of Gases and Liquids, 4th edition,
 * McGraw-Hill, 1987, pp. 43-44
 *
 * and
 *
 * J.E. Killough, et al.: Fifth Comparative Solution Project:
 * Evaluation of Miscible Flood Simulators, Ninth SPE Symposium on
 * Reservoir Simulation, 1987
 */
#ifndef DUMUX_PENG_ROBINSON_PARAMS_MIXTURE_HH
#define DUMUX_PENG_ROBINSON_PARAMS_MIXTURE_HH

#include <algorithm>
#include <dumux/material/constants.hh>

#include "pengrobinsonparams.hh"

namespace Dumux
{

/*!
 * \brief The mixing rule for the oil and the gas phases of the SPE5 problem.
 *
 * This problem comprises \f$H_2O\f$, \f$C_1\f$, \f$C_3\f$, \f$C_6\f$,
 * \f$C_10\f$, \f$C_15\f$ and \f$C_20\f$ as components.
 *
 * See:
 *
 * R. Reid, et al.: The Properties of Gases and Liquids, 4th edition,
 * McGraw-Hill, 1987, pp. 43-44
 *
 * and
 *
 * J.E. Killough, et al.: Fifth Comparative Solution Project:
 * Evaluation of Miscible Flood Simulators, Ninth SPE Symposium on
 * Reservoir Simulation, 1987
 */
template <class Scalar, class FluidSystem, int phaseIdx, bool useSpe5Relations=false>
class PengRobinsonParamsMixture
    : public PengRobinsonParams<Scalar>
{
    enum { numComponents = FluidSystem::numComponents };

    // Peng-Robinson parameters for pure substances
    typedef Dumux::PengRobinsonParams<Scalar> PureParams;

    // the ideal gas constant
    static constexpr Scalar R = Dumux::Constants<Scalar>::R;

public:
    /*!
     * \brief Update Peng-Robinson parameters for the pure components.
     */
    template <class FluidState>
    void updatePure(const FluidState &fluidState)
    {
        updatePure(fluidState.temperature(phaseIdx),
                   fluidState.pressure(phaseIdx));
    }

    /*!
     * \brief Peng-Robinson parameters for the pure components.
     *
     * This method is given by the SPE5 paper.
     */
    void updatePure(Scalar temperature, Scalar pressure)
    {
        Valgrind::CheckDefined(temperature);
        Valgrind::CheckDefined(pressure);

        // Calculate the Peng-Robinson parameters of the pure
        // components
        //
        // See: R. Reid, et al.: The Properties of Gases and Liquids,
        // 4th edition, McGraw-Hill, 1987, p. 43
        for (int i = 0; i < numComponents; ++i) {
            Scalar pc = FluidSystem::criticalPressure(i);
            Scalar omega = FluidSystem::acentricFactor(i);
            Scalar Tr = temperature/FluidSystem::criticalTemperature(i);
            Scalar RTc = R*FluidSystem::criticalTemperature(i);

            Scalar f_omega;

            if (useSpe5Relations) {
                if (omega < 0.49) f_omega = 0.37464  + omega*(1.54226 + omega*(-0.26992));
                else              f_omega = 0.379642 + omega*(1.48503 + omega*(-0.164423 + omega*0.016666));
            }
            else
                f_omega = 0.37464 + omega*(1.54226 - omega*0.26992);

            Valgrind::CheckDefined(f_omega);

            Scalar tmp = 1 + f_omega*(1 - std::sqrt(Tr));
            tmp = tmp*tmp;

            Scalar a = 0.4572355*RTc*RTc/pc * tmp;
            Scalar b = 0.0777961 * RTc / pc;
            assert(std::isfinite(a));
            assert(std::isfinite(b));

            this->pureParams_[i].setA(a);
            this->pureParams_[i].setB(b);
            Valgrind::CheckDefined(this->pureParams_[i].a());
            Valgrind::CheckDefined(this->pureParams_[i].b());
        }

        updateACache_();
    }

    /*!
     * \brief Calculates the "a" and "b" Peng-Robinson parameters for
     *        the mixture.
     *
     * The updatePure() method needs to be called _before_ calling
     * this method!
     */
    template <class FluidState>
    void updateMix(const FluidState &fs)
    {
        Scalar sumx = 0.0;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            sumx += fs.moleFraction(phaseIdx, compIdx);
        sumx = std::max(1e-10, sumx);

        // Calculate the Peng-Robinson parameters of the mixture
        //
        // See: R. Reid, et al.: The Properties of Gases and Liquids,
        // 4th edition, McGraw-Hill, 1987, p. 82
        Scalar a = 0;
        Scalar b = 0;
        for (int compIIdx = 0; compIIdx < numComponents; ++compIIdx) {
            Scalar xi = std::max(0.0, std::min(1.0, fs.moleFraction(phaseIdx, compIIdx)));
            Valgrind::CheckDefined(xi);

            for (int compJIdx = 0; compJIdx < numComponents; ++compJIdx) {
                Scalar xj = std::max(0.0, std::min(1.0, fs.moleFraction(phaseIdx, compJIdx)));
                Valgrind::CheckDefined(xj);

                // mixing rule from Reid, page 82
                a +=  xi * xj * aCache_[compIIdx][compJIdx];

                assert(std::isfinite(a));
            }

            // mixing rule from Reid, page 82
            b += std::max(0.0, xi) * this->pureParams_[compIIdx].b();
            assert(std::isfinite(b));
        }

        assert(b > 0);
        this->setA(a);
        this->setB(b);

        Valgrind::CheckDefined(this->a());
        Valgrind::CheckDefined(this->b());

    }

    /*!
     * \brief Calculates the "a" and "b" Peng-Robinson parameters for
     *        the mixture provided that only a single mole fraction
     *        was changed.
     *
     * The updatePure() method needs to be called _before_ calling
     * this method!
     */
    template <class FluidState>
    void updateSingleMoleFraction(const FluidState &fs,
                                  int compIdx)
    {
        updateMix(fs);
    }

    /*!
     * \brief Return the Peng-Robinson parameters of a pure substance,
     */
    const PureParams &pureParams(int compIdx) const
    { return pureParams_[compIdx]; }

    /*!
     * \brief Returns the Peng-Robinson parameters for a pure component.
     */
    const PureParams &operator[](int compIdx) const
    {
        assert(0 <= compIdx && compIdx < numComponents);
        return pureParams_[compIdx];
    }

    /*!
     * \brief If run under valgrind, this method produces an warning
     *        if the parameters where not determined correctly.
     */
    void checkDefined() const
    {
#ifndef NDEBUG
        for (int i = 0; i < numComponents; ++i)
            pureParams_[i].checkDefined();

        Valgrind::CheckDefined(this->a());
        Valgrind::CheckDefined(this->b());
#endif
    }

protected:
    PureParams pureParams_[numComponents];

private:
    void updateACache_()
    {
        for (int compIIdx = 0; compIIdx < numComponents; ++ compIIdx) {
            for (int compJIdx = 0; compJIdx < numComponents; ++ compJIdx) {
                // interaction coefficient as given in SPE5
                Scalar Psi = FluidSystem::interactionCoefficient(compIIdx, compJIdx);

                aCache_[compIIdx][compJIdx] =
                    std::sqrt(this->pureParams_[compIIdx].a()
                              * this->pureParams_[compJIdx].a())
                    * (1 - Psi);
            }
        }
    }

    Scalar aCache_[numComponents][numComponents];
};

} // end namepace

#endif
