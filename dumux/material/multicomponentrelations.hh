#ifndef MULTICOMPONENTRELATIONS_HH_
#define MULTICOMPONENTRELATIONS_HH_

/**
 * \ingroup material
 * \defgroup properties multicomponent
 */

/** \ingroup properties
 * @brief base class for the computation of multicomponent relations
 */

namespace Dune
{
	
class MultiComp
{
	public:
		/*! \brief solubility of a component (water) in the non-wetting phase
		 *
		 *  \param pN the pressure of the non-wetting phase   
		 *  \param T temperature
		 *  \return the mass fraction (dimensionless)
		 */
		virtual double xWN (double pn, double T=283.15) = 0;
		
		/*! \brief solubility of a component (air) in the wetting phase
		 *
		 *  \param pN the pressure of the non-wetting phase [Pa]   
		 *  \param T temperature [K]
		 *  \return the mass fraction [-]
		 */
		virtual double xAW (double pn, double T=283.15) = 0;
		
		/** @brief Henry coefficient
		 * @param T Temperature \f$ \left[ K \right] \f$
		 * @return Henry coefficient \f$ \left[ 1/Pa \right] \f$
		 */
		virtual double henry (double T=283.15) const = 0;
	
		/*! \brief Antoine equation for calculating the vapor pressure
		 *  \param T temperature [K]
		 *  \return vapor pressure [Pa]
		 */
		virtual double vaporPressure (double T=283.15) const = 0;

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
		double xWN (double pn, double T=283.15)
		{
			double pwsat;
			double result;
			
			pwsat = vaporPressure(T);
			result = pwsat / pn;
				
			result = conversionMoleToMassFraction(result);
				
			return(result);
		}
	  
		double xAW (double pn, double T=283.15)
		{
			double pan;
			double result;
			double hagw;
				
			pan = pn * (1-xWN(pn,T));
			hagw = henry(T);
			result = pan * hagw;
				
			result = conversionMoleToMassFraction(result);
				
			return(result);
	
		}

			  /*! \brief Henry equation for calculating the Henry coefficient
			   *  \param T temperature \f$ \left[ K \right] \f$
			   *  \return Henry Coefficient \f$ \left[ 1/Pa \right] \f$
			   */
		double henry(double T=283.15) const
		{
			double celsius = T - 273.15;
			double result = ((0.8942 + 1.47 * exp(-0.04394*celsius) )*1e-10);
	
			return (result); // [1/Pa]

		}

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

		CWaterAir(const Medium& wP = *(new Uniform), const Medium& nwP = *(new Uniform))
				: MultiComp(wP, nwP)
		{	 }

	protected:
		double conversionMoleToMassFraction(double molefrac) const
		{
			double result;
			double molarMass1 = 0.02896; //air
			double molarMass2 = 0.018016; // water
				
			result = molefrac * molarMass1 / (molarMass1*molefrac + molarMass2*(1-molefrac));
			
			return (result);
		}		
};
}
#endif /*MULTICOMPONENTRELATIONS_HH_*/
