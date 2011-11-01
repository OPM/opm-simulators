/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
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
 * \brief The Peng-Robinson parameters for a mixture
 *
 * See:
 *
 * R. Reid, et al.: The Properties of Gases and Liquids, 4th edition,
 * McGraw-Hill, 1987, pp. 43-44
 */
#ifndef DUMUX_PENG_ROBINSON_PARAMS_MIXTURE_HH
#define DUMUX_PENG_ROBINSON_PARAMS_MIXTURE_HH

#include "pengrobinsonparams.hh"
#include "pengrobinson.hh"

#include <dumux/material/constants.hh>

namespace Dumux
{
/*!
 * \brief The mixing rule for the Peng-Robinson equation of state as given in Reid, p. 82
 *
 * See:
 *
 * R. Reid, et al.: The Properties of Gases and Liquids, 4th edition,
 * McGraw-Hill, 1987, p. 82
 */
template <class Scalar, class StaticParams, int phaseIdx>
class PengRobinsonParamsMixture : public PengRobinsonParams<Scalar>
{
    typedef Dumux::PengRobinsonParams<Scalar> ParentType;

    // Peng-Robinson parameters for pure substances
    typedef Dumux::PengRobinsonParams<Scalar> PureParams;

    // The Peng-Robinson EOS for this mixture
    typedef Dumux::PengRobinson<Scalar> PengRobinson;

    // number of components of which the fluid is composed
    enum { numComponents = StaticParams::numComponents };

    // ideal gas constant
    static constexpr Scalar R = Dumux::Constants<Scalar>::R;

public:
    typedef StaticParams StaticParameters;

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
     * \brief Update Peng-Robinson parameters for the pure components.
     *
     * This method is given by the SPE5 paper.
     */
    void updatePure(Scalar temperature, Scalar pressure)
    {
        // Calculate the Peng-Robinson parameters of the pure
        // components
        //
        // See: R. Reid, et al.: The Properties of Gases and Liquids,
        // 4th edition, McGraw-Hill, 1987, p. 43
        for (int i = 0; i < numComponents; ++i) {
            Scalar pc = StaticParams::criticalPressure(i);
            Scalar omega = StaticParams::acentricFactor(i);
            Scalar Tr = temperature/StaticParams::criticalTemperature(i);
            Scalar RTc = R*StaticParams::criticalTemperature(i);
            Scalar f_omega = 0.37464 + omega*(1.54226 - omega*0.26992);

            Scalar tmp = 1 + f_omega*(1 - std::sqrt(Tr));
            this->pureParams_[i].setA(0.4572355*RTc*RTc/pc
                                      *
                                      tmp*tmp);
            this->pureParams_[i].setB(0.0777961 * RTc / pc);
        }
    }

    /*!
     * \brief Calculates the "a" and "b" Peng-Robinson parameters for
     *        the mixture.
     *
     * The updatePure() method needs to be called _before_ calling
     * this method!
     */
    template <class FluidState>
    void updateMix(const FluidState &fluidState)
    {
        // Calculate the Peng-Robinson parameters of the mixture
        //
        // See: R. Reid, et al.: The Properties of Gases and Liquids,
        // 4th edition, McGraw-Hill, 1987, p. 82
        Scalar a = 0;
        Scalar b = 0;
        for (int i = 0; i < numComponents; ++i) {
            Scalar xi = fluidState.moleFrac(phaseIdx, i);
            for (int j = i; j < numComponents; ++j) {
                Scalar xj = fluidState.moleFrac(phaseIdx, j);

                // interaction coefficient as given in SPE5
                Scalar Psi = StaticParams::interactionCoefficient(i, j);

                // mixing rule from Reid, page 82
                a += xi*xj*std::sqrt(pureParams_[i].a()*pureParams_[j].a())*(1 - Psi);
            }

            // mixing rule from Reid, page 82
            b += xi * pureParams_[i].b();
        }

        this->setA(a);
        this->setB(b);
    }

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
        ParentType::checkDefined();
        for (int i = 0; i < numComponents; ++i)
            pureParams_[i].checkDefined();
#endif
    };

    /*!
     * \brief Return the Peng-Robinson parameters of a pure substance,
     */
    const PureParams &pureParams(int compIdx) const
    { return pureParams_[compIdx]; }


protected:
    PureParams pureParams_[numComponents];
};


} // end namepace

#endif
