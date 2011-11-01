/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
 * \file linearmaterialparams.hh
 *
 * Reference implementation of parameters for the M-phase linear
 * material material.
 */
#ifndef MP_LINEAR_MATERIAL_PARAMS_HH
#define MP_LINEAR_MATERIAL_PARAMS_HH

namespace Dumux
{
/*!
 * \brief Reference implementation of params for the linear M-phase
 *        material material.
 */
template<int numPhasesV, class ScalarT>
class MpLinearMaterialParams
{
public:
    typedef ScalarT Scalar;
    enum { numPhases = numPhasesV };

    MpLinearMaterialParams()
    {
        for (int i = 0; i < numPhases; ++i) {
            setPcMinSat(i, 0.0);
            setPcMaxSat(i, 0.0);
        }
    }

    /*!
     * \brief Return the threshold saturation at which the relative
     *        permeability starts to get regularized.
     *
     * This is simply 10%
     */
    Scalar Sreg(int phaseIdx) const
    { return 0.10; }

    /*!
     * \brief Return the capillary pressure for a phase \f$\alpha\f$ at \f$S_\alpha=0\f$.
     */
    Scalar pcMinSat(int phaseIdx) const
    { return pcMinSat_[phaseIdx]; }

    /*!
     * \brief Set the capillary pressure for a phase \f$\alpha\f$ at \f$S_\alpha=0\f$.
     */
    void setPcMinSat(int phaseIdx, Scalar val)
    { pcMinSat_[phaseIdx] = val; }

    /*!
     * \brief Return the capillary pressure for a phase \f$\alpha\f$ at \f$S_\alpha=1\f$.
     */
    Scalar pcMaxSat(int phaseIdx) const
    { return pcMaxSat_[phaseIdx]; }

    /*!
     * \brief Set the capillary pressure for a phase \f$\alpha\f$ at \f$S_\alpha=1\f$.
     */
    void setPcMaxSat(int phaseIdx, Scalar val)
    { pcMaxSat_[phaseIdx] = val; }

    /*!
     * \brief Return the threshold saturation respective phase below
     *        which the relative permeability gets regularized.
     *
     * This is just 5%. If you need a different value, write your own
     * parameter class.
     */
    Scalar krLowS(int phaseIdx) const
    { return 0.05; }

    /*!
     * \brief Return the threshold saturation of the respective phase
     *        above which the relative permeability gets regularized.
     *
     * This is just 95%. If you need a different value, write your own
     * parameter class.
     */
    Scalar krHighS(int phaseIdx) const
    { return 0.95; }

private:
    Scalar pcMaxSat_[numPhases];
    Scalar pcMinSat_[numPhases];
};
} // namespace Dumux

#endif
