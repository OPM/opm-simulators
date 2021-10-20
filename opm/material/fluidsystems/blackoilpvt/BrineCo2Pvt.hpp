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
 * \copydoc Opm::BrineCo2Pvt
 */
#ifndef OPM_BRINE_CO2_PVT_HPP
#define OPM_BRINE_CO2_PVT_HPP

#include <opm/material/Constants.hpp>

#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/components/Brine.hpp>
#include <opm/material/components/SimpleHuDuanH2O.hpp>
#include <opm/material/components/CO2.hpp>
#include <opm/material/common/UniformTabulated2DFunction.hpp>
#include <opm/material/components/TabulatedComponent.hpp>
#include <opm/material/binarycoefficients/H2O_CO2.hpp>
#include <opm/material/binarycoefficients/Brine_CO2.hpp>
#include <opm/material/components/co2tables.inc>


#if HAVE_ECL_INPUT
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#endif

#include <vector>

namespace Opm {
/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the liquid phase
 * for a CO2-Brine system
 */
template <class Scalar>
class BrineCo2Pvt
{
    typedef std::vector<std::pair<Scalar, Scalar> > SamplingPoints;
    static const bool extrapolate = true;
    //typedef H2O<Scalar> H2O_IAPWS;
    //typedef Brine<Scalar, H2O_IAPWS> Brine_IAPWS;
    //typedef TabulatedComponent<Scalar, H2O_IAPWS> H2O_Tabulated;
    //typedef TabulatedComponent<Scalar, Brine_IAPWS> Brine_Tabulated;

    //typedef H2O_Tabulated H2O;
    //typedef Brine_Tabulated Brine;


public:
    typedef SimpleHuDuanH2O<Scalar> H2O;
    typedef ::Opm::Brine<Scalar, H2O> Brine;
    typedef ::Opm::CO2<Scalar, CO2Tables> CO2;

    typedef Tabulated1DFunction<Scalar> TabulatedOneDFunction;

    //! The binary coefficients for brine and CO2 used by this fluid system
    typedef BinaryCoeff::Brine_CO2<Scalar, H2O, CO2> BinaryCoeffBrineCO2;

    explicit BrineCo2Pvt() = default;
    BrineCo2Pvt(const std::vector<Scalar>& brineReferenceDensity,
                const std::vector<Scalar>& co2ReferenceDensity,
                const std::vector<Scalar>& salinity)
        : brineReferenceDensity_(brineReferenceDensity),
          co2ReferenceDensity_(co2ReferenceDensity),
          salinity_(salinity)
    {
    }
#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize the parameters for Brine-CO2 system using an ECL deck.
     *
     */
    void initFromState(const EclipseState& eclState, const Schedule&)
    {
        if( !eclState.getTableManager().getDensityTable().empty()) {
            std::cerr << "WARNING: CO2STOR is enabled but DENSITY is in the deck. \n" <<
                         "The surface density is computed based on CO2-BRINE PVT at standard conditions (STCOND) and DENSITY is ignored " << std::endl;
        }

        if( eclState.getTableManager().hasTables("PVDO") || !eclState.getTableManager().getPvtoTables().empty()) {
            std::cerr << "WARNING: CO2STOR is enabled but PVDO or PVTO is in the deck. \n" <<
                         "BRINE PVT properties are computed based on the Hu et al. pvt model and PVDO/PVTO input is ignored. " << std::endl;
        }

        // We only supported single pvt region for the co2-brine module
        size_t numRegions = 1;
        setNumRegions(numRegions);
        size_t regionIdx = 0;
        // Currently we only support constant salinity
        const Scalar molality = eclState.getTableManager().salinity(); // mol/kg
        const Scalar MmNaCl = 58e-3; // molar mass of NaCl [kg/mol]
        // convert to mass fraction
        Brine::salinity = 1 / ( 1 + 1 / (molality*MmNaCl)); //
        salinity_[regionIdx] = Brine::salinity;
        // set the surface conditions using the STCOND keyword
        Scalar T_ref = eclState.getTableManager().stCond().temperature;
        Scalar P_ref = eclState.getTableManager().stCond().pressure;

        brineReferenceDensity_[regionIdx] = Brine::liquidDensity(T_ref, P_ref, extrapolate);
        co2ReferenceDensity_[regionIdx] = CO2::gasDensity(T_ref, P_ref, extrapolate);
    }
#endif

    void setNumRegions(size_t numRegions)
    {
        brineReferenceDensity_.resize(numRegions);
        co2ReferenceDensity_.resize(numRegions);
        salinity_.resize(numRegions);
    }


    /*!
     * \brief Initialize the reference densities of all fluids for a given PVT region
     */
    void setReferenceDensities(unsigned regionIdx,
                               Scalar rhoRefBrine,
                               Scalar rhoRefCO2,
                               Scalar /*rhoRefWater*/)
    {
        brineReferenceDensity_[regionIdx] = rhoRefBrine;
        co2ReferenceDensity_[regionIdx] = rhoRefCO2;
    }


    /*!
     * \brief Finish initializing the oil phase PVT properties.
     */
    void initEnd()
    {

    }

    /*!
     * \brief Return the number of PVT regions which are considered by this PVT-object.
     */
    unsigned numRegions() const
    { return brineReferenceDensity_.size(); }

    /*!
     * \brief Returns the specific enthalpy [J/kg] of gas given a set of parameters.
     */
    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx,
                        const Evaluation& temperature,
                        const Evaluation& pressure,
                        const Evaluation& Rs) const
    {

        const Evaluation xlCO2 = convertXoGToxoG_(convertRsToXoG_(Rs,regionIdx));
        return (liquidEnthalpyBrineCO2_(temperature,
                                       pressure,
                                       salinity_[regionIdx],
                                       xlCO2)
        - pressure / density_(regionIdx, temperature, pressure, Rs));
    }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& /*Rs*/) const
    {
        //TODO: The viscosity does not yet depend on the composition
        return saturatedViscosity(regionIdx, temperature, pressure);
    }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of oil saturated gas at given pressure.
     */
    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned /*regionIdx*/,
                                  const Evaluation& temperature,
                                  const Evaluation& pressure) const
    {
        return Brine::liquidViscosity(temperature, pressure);
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& temperature,
                                            const Evaluation& pressure,
                                            const Evaluation& Rs) const
    {
        return density_(regionIdx, temperature, pressure, Rs)/brineReferenceDensity_[regionIdx];
    }

    /*!
     * \brief Returns the formation volume factor [-] of brine saturated with CO2 at a given pressure.
     */
    template <class Evaluation>
    Evaluation saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                                     const Evaluation& temperature,
                                                     const Evaluation& pressure) const
    {
        Evaluation rsSat = rsSat_(regionIdx, temperature, pressure);
        return density_(regionIdx, temperature, pressure, rsSat)/brineReferenceDensity_[regionIdx];
    }

    /*!
     * \brief Returns the saturation pressure of the brine phase [Pa]
     *        depending on its mass fraction of the gas component
     *
     * \param Rs
     */
    template <class Evaluation>
    Evaluation saturationPressure(unsigned /*regionIdx*/,
                                  const Evaluation& /*temperature*/,
                                  const Evaluation& /*Rs*/) const
    {
        throw std::runtime_error("Requested the saturation pressure for the brine-co2 pvt module. Not yet implemented.");
    }

    /*!
     * \brief Returns the gas dissoluiton factor \f$R_s\f$ [m^3/m^3] of the liquid phase.
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure,
                                             const Evaluation& /*oilSaturation*/,
                                             const Evaluation& /*maxOilSaturation*/) const
    {
        //TODO support VAPPARS
        return rsSat_(regionIdx, temperature, pressure);
    }

    /*!
     * \brief Returns thegas dissoluiton factor  \f$R_s\f$ [m^3/m^3] of the liquid phase.
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure) const
    {
        return rsSat_(regionIdx, temperature, pressure);
    }

    const Scalar oilReferenceDensity(unsigned regionIdx) const
    { return brineReferenceDensity_[regionIdx]; }

    const Scalar gasReferenceDensity(unsigned regionIdx) const
    { return co2ReferenceDensity_[regionIdx]; }

    const Scalar salinity(unsigned regionIdx) const
    { return salinity_[regionIdx]; }

    bool operator==(const BrineCo2Pvt<Scalar>& data) const
    {
        return co2ReferenceDensity_ == data.co2ReferenceDensity_ &&
                brineReferenceDensity_ == data.brineReferenceDensity_;
    }

    template <class Evaluation>
    Evaluation diffusionCoefficient(const Evaluation& temperature,
                                    const Evaluation& pressure,
                                    unsigned /*compIdx*/) const
    {
        //Diffusion coefficient of CO2 in pure water according to (McLachlan and Danckwerts, 1972)
        const Evaluation log_D_H20 = -4.1764 + 712.52 / temperature - 2.5907e5 / (temperature*temperature);

        //Diffusion coefficient of CO2 in the brine phase modified following (Ratcliff and Holdcroft,1963 and Al-Rawajfeh, 2004)
        const Evaluation& mu_H20 = H2O::liquidViscosity(temperature, pressure, extrapolate); // Water viscosity
        const Evaluation& mu_Brine = Brine::liquidViscosity(temperature, pressure); // Brine viscosity
        const Evaluation log_D_Brine = log_D_H20 - 0.87*log10(mu_Brine / mu_H20);

        return pow(Evaluation(10), log_D_Brine) * 1e-4; // convert from cm2/s to m2/s
    }

private:
    std::vector<Scalar> brineReferenceDensity_;
    std::vector<Scalar> co2ReferenceDensity_;
    std::vector<Scalar> salinity_;

    template <class LhsEval>
    LhsEval density_(unsigned regionIdx,
                     const LhsEval& temperature,
                     const LhsEval& pressure,
                     const LhsEval& Rs) const
    {
        LhsEval xlCO2 = convertXoGToxoG_(convertRsToXoG_(Rs,regionIdx));
        LhsEval result = liquidDensity_(temperature,
                                        pressure,
                                        xlCO2);

        Valgrind::CheckDefined(result);
        return result;
    }


    template <class LhsEval>
    LhsEval liquidDensity_(const LhsEval& T,
                           const LhsEval& pl,
                           const LhsEval& xlCO2) const
    {
        Valgrind::CheckDefined(T);
        Valgrind::CheckDefined(pl);
        Valgrind::CheckDefined(xlCO2);

        if(!extrapolate && T < 273.15) {
            std::ostringstream oss;
            oss << "Liquid density for Brine and CO2 is only "
                   "defined above 273.15K (is "<<T<<"K)";
            throw NumericalIssue(oss.str());
        }
        if(!extrapolate && pl >= 2.5e8) {
            std::ostringstream oss;
            oss << "Liquid density for Brine and CO2 is only "
                   "defined below 250MPa (is "<<pl<<"Pa)";
            throw NumericalIssue(oss.str());
        }

        const LhsEval& rho_brine = Brine::liquidDensity(T, pl, extrapolate);
        const LhsEval& rho_pure = H2O::liquidDensity(T, pl, extrapolate);
        const LhsEval& rho_lCO2 = liquidDensityWaterCO2_(T, pl, xlCO2);
        const LhsEval& contribCO2 = rho_lCO2 - rho_pure;

        return rho_brine + contribCO2;
    }

    template <class LhsEval>
    LhsEval liquidDensityWaterCO2_(const LhsEval& temperature,
                                          const LhsEval& pl,
                                          const LhsEval& xlCO2) const
    {
        Scalar M_CO2 = CO2::molarMass();
        Scalar M_H2O = H2O::molarMass();

        const LhsEval& tempC = temperature - 273.15;        /* tempC : temperature in °C */
        const LhsEval& rho_pure = H2O::liquidDensity(temperature, pl, extrapolate);
        // calculate the mole fraction of CO2 in the liquid. note that xlH2O is available
        // as a function parameter, but in the case of a pure gas phase the value of M_T
        // for the virtual liquid phase can become very large
        const LhsEval xlH2O = 1.0 - xlCO2;
        const LhsEval& M_T = M_H2O * xlH2O + M_CO2 * xlCO2;
        const LhsEval& V_phi =
            (37.51 +
             tempC*(-9.585e-2 +
                    tempC*(8.74e-4 -
                           tempC*5.044e-7))) / 1.0e6;
        return 1/ (xlCO2 * V_phi/M_T + M_H2O * xlH2O / (rho_pure * M_T));
    }

    /*!
     * \brief Convert a gas dissolution factor to the the corresponding mass fraction
     *        of the gas component in the oil phase.
     */
    template <class LhsEval>
    LhsEval convertRsToXoG_(const LhsEval& Rs, unsigned regionIdx) const
    {
        Scalar rho_oRef = brineReferenceDensity_[regionIdx];
        Scalar rho_gRef = co2ReferenceDensity_[regionIdx];

        const LhsEval& rho_oG = Rs*rho_gRef;
        return rho_oG/(rho_oRef + rho_oG);
    }


    /*!
     * \brief Convert a gas mass fraction in the oil phase the corresponding mole fraction.
     */
    template <class LhsEval>
    LhsEval convertXoGToxoG_(const LhsEval& XoG) const
    {
        Scalar M_CO2 = CO2::molarMass();
        Scalar M_Brine = Brine::molarMass();
        return XoG*M_Brine / (M_CO2*(1 - XoG) + XoG*M_Brine);
    }


    /*!
     * \brief Convert a gas mole fraction in the oil phase the corresponding mass fraction.
     */
    template <class LhsEval>
    LhsEval convertxoGToXoG(const LhsEval& xoG) const
    {
        Scalar M_CO2 = CO2::molarMass();
        Scalar M_Brine = Brine::molarMass();

        return xoG*M_CO2 / (xoG*(M_CO2 - M_Brine) + M_Brine);
    }


    /*!
     * \brief Convert the mass fraction of the gas component in the oil phase to the
     *        corresponding gas dissolution factor.
     */
    template <class LhsEval>
    LhsEval convertXoGToRs(const LhsEval& XoG, unsigned regionIdx) const
    {
        Scalar rho_oRef = brineReferenceDensity_[regionIdx];
        Scalar rho_gRef = co2ReferenceDensity_[regionIdx];

        return XoG/(1.0 - XoG)*(rho_oRef/rho_gRef);
    }


    template <class LhsEval>
    LhsEval rsSat_(unsigned regionIdx,
                   const LhsEval& temperature,
                   const LhsEval& pressure) const
    {
        // calulate the equilibrium composition for the given
        // temperature and pressure. 
        LhsEval xgH2O;
        LhsEval xlCO2;
        BinaryCoeffBrineCO2::calculateMoleFractions(temperature,
                                                    pressure,
                                                    salinity_[regionIdx],
                                                    /*knownPhaseIdx=*/-1,
                                                    xlCO2,
                                                    xgH2O,
                                                    extrapolate);

        // normalize the phase compositions
        xlCO2 = max(0.0, min(1.0, xlCO2));

        return convertXoGToRs(convertxoGToXoG(xlCO2), regionIdx);
    }

    template <class LhsEval>
    static LhsEval liquidEnthalpyBrineCO2_(const LhsEval& T,
                                           const LhsEval& p,
                                           Scalar S, // salinity
                                           const LhsEval& X_CO2_w)
    {
        /* X_CO2_w : mass fraction of CO2 in brine */

        /* same function as enthalpy_brine, only extended by CO2 content */

        /*Numerical coefficents from PALLISER*/
        static Scalar f[] = {
            2.63500E-1, 7.48368E-6, 1.44611E-6, -3.80860E-10
        };

        /*Numerical coefficents from MICHAELIDES for the enthalpy of brine*/
        static Scalar a[4][3] = {
            { 9633.6, -4080.0, +286.49 },
            { +166.58, +68.577, -4.6856 },
            { -0.90963, -0.36524, +0.249667E-1 },
            { +0.17965E-2, +0.71924E-3, -0.4900E-4 }
        };

        LhsEval theta, h_NaCl;
        LhsEval h_ls1, d_h;
        LhsEval delta_h;
        LhsEval delta_hCO2, hg, hw;

        theta = T - 273.15;

        // Regularization
        Scalar scalarTheta = scalarValue(theta);
        Scalar S_lSAT = f[0] + scalarTheta*(f[1] + scalarTheta*(f[2] + scalarTheta*f[3]));
        if (S > S_lSAT)
            S = S_lSAT;

        hw = H2O::liquidEnthalpy(T, p) /1E3; /* kJ/kg */

        /*DAUBERT and DANNER*/
        /*U=*/h_NaCl = (3.6710E4*T + 0.5*(6.2770E1)*T*T - ((6.6670E-2)/3)*T*T*T
                        +((2.8000E-5)/4)*(T*T*T*T))/(58.44E3)- 2.045698e+02; /* kJ/kg */

        Scalar m = 1E3/58.44 * S/(1-S);
        int i = 0;
        int j = 0;
        d_h = 0;

        for (i = 0; i<=3; i++) {
            for (j=0; j<=2; j++) {
                d_h = d_h + a[i][j] * pow(theta, static_cast<Scalar>(i)) * std::pow(m, j);
            }
        }
        /* heat of dissolution for halite according to Michaelides 1971 */
        delta_h = (4.184/(1E3 + (58.44 * m)))*d_h;

        /* Enthalpy of brine without CO2 */
        h_ls1 =(1-S)*hw + S*h_NaCl + S*delta_h; /* kJ/kg */

        /* heat of dissolution for CO2 according to Fig. 6 in Duan and Sun 2003. (kJ/kg)
           In the relevant temperature ranges CO2 dissolution is
           exothermal */
        delta_hCO2 = (-57.4375 + T * 0.1325) * 1000/44;

        /* enthalpy contribution of CO2 (kJ/kg) */
        hg = CO2::gasEnthalpy(T, p, extrapolate)/1E3 + delta_hCO2;

        /* Enthalpy of brine with dissolved CO2 */
        return (h_ls1 - X_CO2_w*hw + hg*X_CO2_w)*1E3; /*J/kg*/
    }

};

} // namespace Opm

#endif
