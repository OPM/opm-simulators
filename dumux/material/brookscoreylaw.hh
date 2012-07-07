#ifndef DUNE_BROOKSCOREYLAW_HH
#define DUNE_BROOKSCOREYLAW_HH

#include "dumux/material/twophaserelations.hh"

namespace Dune
{
	/*!\ingroup material 
	 * \brief Brooks-Corey mobility/saturation relation.  
	 *
	 *  Employs the Brooks-Corey non-linear relative permeability/saturation relation, namely, 
	 *  \f{align*} 
	 *  S_\text{e} = \frac{S - S_\text{r}}{1 - S_\text{r}}, \quad 
	 *  \lambda_\text{w} = \mu_\text{w}^{-1} S_\text{e}^{\frac{2 + 3\lambda}{\lambda}},  \quad 
	 *  \lambda_\text{n} = \mu_\text{n}^{-1} \left( 1 - S_\text{e}^{\frac{2 + 2\lambda}{\lambda}}\right)
	 *  \left( 1 - S_\text{e} \right)^2, \quad
	 *  \lambda = \lambda_\text{w} + \lambda_\text{n}. 
	 *  \f} 
	 */
	class BrooksCoreyLaw : public TwoPhaseRelations
	{
	public:
	  double pC (double saturationW) 
	  {
		  double Se = (saturationW - this->wettingPhase.Sr())
						/(1. - this->wettingPhase.Sr() - this->nonwettingPhase.Sr());

		  if (Se > epsPC) 
			  return (p0*pow(Se, -1.0/lambda));
		  else {
			  double dpCEps = dPdS(epsPC,283.15, 1e5);
			  return (dpCEps*(Se - epsPC) + p0*pow(epsPC, -1.0/lambda));
		  }
	  }
	
	  double pC (double saturationW, const FieldVector<double, 4>& parameters) 
	  {
		  DUNE_THROW(NotImplemented,"BrooksCoreyLaw :: pC (double, const FieldVector<double, 4>&");
	  }
	  
	  double dPdS (double saturationW, double T=283.15, double p=1e5)
	  {
		  double Se = (saturationW - this->wettingPhase.Sr())
			/(1. - this->wettingPhase.Sr() - this->nonwettingPhase.Sr());

		  return (p0/lambda*pow(std::max(Se, epsPC), -(1.0 + lambda)/lambda));
	  }
	  
	  double saturationW (double pC, double T=283.15) 
	  {
		  return ((1 - this->wettingPhase.Sr())*pow(p0/pC, lambda) + this->wettingPhase.Sr());
	  }

	  double dSdP (double pC, double T=283.15) const 
	  {
		  return (-(1 - this->wettingPhase.Sr())*pow(p0/pC, lambda-1)*lambda*pow(1.0/pC, 2));
	  }

	  BrooksCoreyLaw(const Medium& wP = *(new Uniform), const Medium& nwP = *(new Uniform), 
			  double l = 0.0, double p = 0.0, double eps = 0.0) 
	    : TwoPhaseRelations(wP, nwP, false)
	  {
	    lambda = l ? l : 2.0;
	    p0 = p ? p : 1.0e5;
	    epsPC = eps ? eps : 0.01;
	  }
	  
	protected:
		double krw (double saturationW) const
		{
		    double Se = (saturationW - this->wettingPhase.Sr())
		    			/(1. - this->wettingPhase.Sr() - this->nonwettingPhase.Sr());
		    double exponent = (2. + 3*lambda) / lambda;
		    return pow(Se, exponent);
		}
	  
		double krn (double saturationN) const 
		{
		    double Se = ((1.-saturationN) - this->wettingPhase.Sr())
		    			/(1. - this->wettingPhase.Sr() - this->nonwettingPhase.Sr());
		    double exponent = (2. + lambda) / lambda;
		    return pow(1.-Se, 2) * ( 1. - pow(Se, exponent) );
		}
		
		double lambda; //!< material parameter for the relative permeability/saturation relation
		double p0; //!< capillary entry pressure
		double epsPC; //!< treshold for linearization of capillary pressure
	};
}

#endif


