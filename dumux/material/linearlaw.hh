#ifndef DUNE_LINEARLAW_HH
#define DUNE_LINEARLAW_HH

#include "dumux/material/twophaserelations.hh"

namespace Dune
{
	/** \ingroup material
	 * @brief Represents the linear relative permability - saturation relation \f$ k_{r\alpha}(S) = \frac{S-S_{r\alpha}}{1-S_{r\alpha}} \f$  
	 */
	class LinearLaw : public TwoPhaseRelations 
	{
	public:
	  double pC (double saturationW) 
	  {
		  return 0.0;
	  }
	
		double pC (double saturationW, const FieldVector<double, 4>& parameters) 
		{
			return 0;
		}

		double dPdS (double saturationW, double T=283.15, double p=1e5) 
	  {
		  return 0.0;  
	  }
	  
	  double saturationW (double pC, double T=283.15) 
	  {
		  return (0.0);
	  }

	  double dSdP (double pC, double T=283.15) const 
	  {
		  return (0.0);
	  }

	  LinearLaw(const Medium& wP = *(new Uniform), const Medium& nwP = *(new Uniform))
	  : TwoPhaseRelations(wP, nwP, true)
	  {	 }
	  
	protected:
		double krw (double saturationW) const
		{
		    return std::max(std::min(saturationW - this->wettingPhase.Sr(), 1.), 0.);			
		}
	  
		double krw (double saturationW, const FieldVector<double, 4>& parameters) const 
		{
			return krw(saturationW);
		}

		double krn (double saturationN) const 
		{
		    return std::max( std::min(saturationN - this->nonwettingPhase.Sr(), 1.), 0.);
		}

		double krn (double saturationN, const FieldVector<double, 4>& parameters) const 
		{
			return krn(saturationN);
		}
	};
}

#endif


