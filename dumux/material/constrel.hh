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
		
		virtual double Xwn (double pn=1e5, double T=283.15)
		{
			double pwsat;
			double result;
		
			pwsat = vaporPressure(T);
			result = pwsat / pn;
			
			return(result);
		}
	
		virtual double Xaw (double pn=1e5, double T=283.15)
		{
			double pan;
			double result;
			double hagw;
			
			pan = pn * (1-Xwn(pn,T));
			hagw = HenryWaterAir(T);
			result = pan * hagw;
			
			return(result);
		}

		virtual ~Solubility()
		{}
		
		
	protected:

		virtual double HenryWaterAir(double T=283.15)
		{
			double celsius = T - 273.15;
			double result = ((0.8942 + 1.47 * exp(-0.04394*celsius) )*1e-10);
	
			return (result); // [1/Pa]
		}
		
		virtual double vaporPressure(double T=283.15)//Antoine(double T=283.15)
		{
			const double constA = 8.19621; 
			const double constB = 1730.63;
			const double constC = 233.436;
			double celsius;
			double exponent, pwsat;
	
			celsius = T - 273.15;
	
			exponent = constA - (constB / (celsius + constC));
	
			pwsat = pow (10.0, exponent) *100; /*1mbar = 100Pa*/
	
			return(pwsat);
		}

};

}
#endif

