#ifndef DUNE_CONSTREL_HH
#define DUNE_CONSTREL_HH
/**
 * \ingroup material
 * \defgroup properties Constitutive relationships
 */

/** \ingroup properties
 * @brief base class for constitutive relationships
 */
namespace Dune
{
	
class Solubility
{
	public:
	//	virtual double operator() (double T=283.15) = 0; 
		
		//mole fractions, converted in mass fractions
		virtual double Xwn (double pn=1e5, double T=283.15)
		{
			double pwsat;
			double result;
		
			pwsat = vaporPressure(T);
			result = pwsat / pn;
			
			result = conversionMoleToMassFraction(result);
			
			return(result);
		}
	
		virtual double Xaw (double pn=1e5, double T=283.15)
		{
			double pan;
			double result;
			double hagw;
			
			pan = pn * (1-Xwn(pn,T));
			hagw = henryWaterAir(T);
			result = pan * hagw;
			
			result = conversionMoleToMassFraction(result);
			
			return(result);
		}

		virtual ~Solubility()
		{}
		
		
	protected:

		virtual double henryWaterAir(double T=283.15)
		{
			double celsius = T - 273.15;
			double result = ((0.8942 + 1.47 * exp(-0.04394*celsius) )*1e-10);
	
			return (result); // [1/Pa]
		}

		// Antoine equation for calculating the vapor pressure
		virtual double vaporPressure(double T=283.15)
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

		double conversionMoleToMassFraction(double molefrac)
		{
			double result;
			double molarMass1 = 0.02896; //air
			double molarMass2 = 0.018016; // water
				
			result = molefrac * molarMass1 / (molarMass1*molefrac + molarMass2*(1-molefrac));
			
			return (result);
		}
};

}
#endif

