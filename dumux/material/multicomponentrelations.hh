#ifndef MULTICOMPONENTRELATIONS_HH_
#define MULTICOMPONENTRELATIONS_HH_

/**
 * \ingroup material
 * \defgroup properties multicomponent
 * \author Klaus Mosthaf
 */


namespace Dune
{

/** \ingroup properties
 *  \brief base class for the computation of multicomponent relations
 */
class MultiComp
{
	public:
		/*! \brief solubility of a component (water) in the non-wetting phase
		 *	\param pn non-wetting phase pressure \f$ \left[ Pa \right] \f$ 
		 *  \param T temperature
		 *  \return the mass fraction of water in the non-wetting phase \f$ \left[ kg/kg \right] \f$
		 */
		virtual double xWN (const double pn, double T=283.15) = 0;
		
		/*! \brief solubility of a component (air) in the wetting phase
		 *	\param pn non-wetting phase pressure \f$ \left[ Pa \right] \f$ 
		 *  \param T temperature [K]
		 *  \return mass fraction of gas in the wetting phase \f$ \left[ kg/kg \right] \f$
		 */
		virtual double xAW (const double pn, double T=283.15) = 0;
		
		/*! \brief solubility of a component (water) in the non-wetting phase
		 *	\param pn non-wetting phase pressure \f$ \left[ Pa \right] \f$ 
		 *  \param T temperature
		 *  \return the mass fraction of water in the non-wetting phase \f$ \left[ mol/mol \right] \f$
		 */
		virtual double xWNmolar (const double pn, double T=283.15) = 0;
		
		/*! \brief solubility of a component (air) in the wetting phase
		 *	\param pn non-wetting phase pressure \f$ \left[ Pa \right] \f$ 
		 *  \param T temperature [K]
		 *  \return the mass fraction of water in the non-wetting phase \f$ \left[ mol/mol \right] \f$
		 */
		virtual double xAWmolar (const double pn, double T=283.15) = 0;
		
		/** @brief Henry coefficient
		 * \param pn non-wetting phase pressure \f$ \left[ Pa \right] \f$ 
		 * @param T Temperature \f$ \left[ K \right] \f$
		 * @return Henry coefficient \f$ \left[ 1/Pa \right] \f$
		 */
		virtual double henry (double T=283.15) const = 0;
	
		/*! \brief Antoine equation for calculating the vapor pressure
		 *  \param T temperature [K]
		 *  \return vapor pressure [Pa]
		 */
		virtual double vaporPressure (double T=283.15) const = 0;

		/*! \brief converts mole fractions into mass fractions
		 *  \param massfrac mole fraction [mol/mol]
		 *  \param phase phase for which a conversion is to be done [-]
		 *  \return mass fraction [kg/kg]
		 */
		virtual double conversionMoleToMassFraction(double massfrac, int phase) const = 0;

		/*! \brief converts mass fractions into mole fractions
		 *  \param massfrac mass fraction [kg/kg]
		 *  \param phase phase for which a conversion is to be done [-]
		 *  \return mole fraction [mol/mol]
		 */
		virtual double conversionMassToMoleFraction(double massfrac, int phase) const = 0;


		MultiComp(const Medium& wP = *(new Uniform), const Medium& nwP = *(new Uniform))
				: wettingPhase(wP), nonwettingPhase(nwP)
		{	 }

		virtual ~MultiComp()
		{}

		const Medium& wettingPhase; //!< contains the properties of the wetting phase 
		const Medium& nonwettingPhase; //!< contains the properties of the nonwetting phase 

};
		
		
class CWaterAir : public MultiComp
{
	public:
		/*! \brief equation for calculating the mass fraction in the nonwetting phase
		 *	\param pn non-wetting phase pressure \f$ \left[ Pa \right] \f$ 
		 *  \param T temperature \f$ \left[ K \right] \f$
		 *  \return mass fraction \f$ \left[ kg/kg \right] \f$
		 */
		double xWN (const double pn, const double T=283.15)
		{
			double pwsat;
			double result;
			
			pwsat = vaporPressure(T);
			result = pwsat / pn;
				
			result = conversionMoleToMassFraction(result, 1);
				
			return(result);
		}
	  

		/*! \brief equation for calculating the mass fraction in the wetting phase
		 *	\param pn non-wetting phase pressure \f$ \left[ Pa \right] \f$ 
		 *  \param T temperature \f$ \left[ K \right] \f$
		 *  \return mass fraction \f$ \left[ kg/kg \right] \f$
		 */
		double xAW (const double pn, const double T=283.15)
		{
			double pan;
			double result;
			double hagw;
				
			pan = pn * (1-xWNmolar(pn,T));
			hagw = henry(T);
			result = pan * hagw;
				
			result = conversionMoleToMassFraction(result, 0);
				
			return(result);
		}
		
		/*! \brief equation for calculating the mole fraction in the nonwetting phase
		 *	\param pn non-wetting phase pressure \f$ \left[ Pa \right] \f$ 
		 *  \param T temperature \f$ \left[ K \right] \f$
		 *  \return mole fraction \f$ \left[ mol/mol \right] \f$
		 */
		double xWNmolar (const double pn, const double T=283.15)
		{
			double pwsat;
			double result;
			
			pwsat = vaporPressure(T);
			result = pwsat / pn;
				
			return(result);
		}
		
		/*! \brief equation for calculating the mole fraction in the wetting phase
		 *	\param pn non-wetting phase pressure \f$ \left[ Pa \right] \f$ 
		 *  \param T temperature \f$ \left[ K \right] \f$
		 *  \return mole fraction \f$ \left[ mol/mol \right] \f$
		 */
		double xAWmolar (const double pn, const double T=283.15)
		{
			double pan;
			double result;
			double hagw;
				
			pan = pn * (1-xWNmolar(pn,T));
			hagw = henry(T);
			result = pan * hagw;
				
			return(result);
		}

		/*! \brief equation for calculating the inverse Henry coefficient
		 *  \param T temperature \f$ \left[ K \right] \f$
		 *  \return Henry Coefficient \f$ \left[ 1/Pa \right] \f$
		 */
		double henry(double T=283.15) const
		{
			double celsius = T - 273.15;
			double result = ((0.8942 + 1.47 * exp(-0.04394*celsius) )*1e-10);
	
			return (result); // [1/Pa]
		}

		/** @brief calculates vapor pressure
		 *  @param T temperature \f$ \left[ K \right] \f$
		 *  @return vapor pressure \f$ \left[ Pa \right] \f$
		 */
		double vaporPressure(double T=283.15) const
		{
			const double constA = 8.19621; 
			const double constB = 1730.63;
			const double constC = 233.426;

			double celsius;
			double exponent, psat;
	
			celsius = T - 273.15;
	
			exponent = constA - (constB / (celsius + constC));
	
			psat = pow (10.0, exponent) *100; //1mbar = 100Pa
	
			return(psat);
		}

		/** @brief converts mole fractions into mass fractions
		 */
		double conversionMoleToMassFraction(double molefrac, int phase) const
		{
			enum {wPhase = 0, nPhase = 1};
			
			double result;
			double molarMass1, molarMass2;

			if (phase == wPhase){			
				molarMass1 = this->wettingPhase.molarMass();
				molarMass2 = this->nonwettingPhase.molarMass();
			}
			if (phase == nPhase){			
				molarMass1 = this->nonwettingPhase.molarMass();
				molarMass2 = this->wettingPhase.molarMass();
			}

			result = molefrac * molarMass1 / (molarMass1*molefrac + molarMass2*(1-molefrac));
			
			return (result);
		}		

		/** @brief converts mass fractions into mole fractions
		 */
		double conversionMassToMoleFraction(double massfrac, int phase) const
		{
			enum {wPhase = 0, nPhase = 1};
			
			double result;
			double molarMass1, molarMass2;

			if (phase == wPhase){			
				molarMass1 = this->wettingPhase.molarMass();
				molarMass2 = this->nonwettingPhase.molarMass();
			}
			if (phase == nPhase){			
				molarMass1 = this->nonwettingPhase.molarMass();
				molarMass2 = this->wettingPhase.molarMass();
			}

			result = massfrac * molarMass2 / (molarMass1*(1-massfrac) + molarMass2*massfrac);
			
			return (result);
		}		

		CWaterAir(const Medium& wP = *(new Uniform), const Medium& nwP = *(new Uniform))
				: MultiComp(wP, nwP)
		{	 }

		
};
}
#endif /*MULTICOMPONENTRELATIONS_HH_*/
