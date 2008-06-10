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
      double Swr = parameters[0];
      double Snr = parameters[1];
      lambda = parameters[2];
      p0 = parameters[3];
      double maxpc = 1e20;
         
      double Se,dSedSw,slope,offset;

      //effective Saturation
      Se = (saturationW - Swr)/(1-Snr-Swr);
      if(Se<0.0) Se = 0.0;
      if(Se>1.0) Se = 1.0;

      if (Se > epsPC) 
	return (std::min(p0*pow(Se, -1.0/lambda),maxpc));
      else {
	dSedSw=1/(1-Snr-Swr);
	slope=-p0*dSedSw/pow(epsPC,(1+1/lambda))/lambda;
	offset=p0*pow(epsPC,(-1/lambda));
	return (std::min((Se - epsPC)*slope +offset,maxpc));
      }
    }
	  
    double dPdS (double saturationW, double T=283.15, double p=1e5)
    {
      double Swr = this->wettingPhase.Sr();
      double Snr = this->nonwettingPhase.Sr();
      
      double Se = (saturationW - Swr)
	/(1. - Swr - Snr);
      
	if(Se<0.0) Se = 0.0;
	if(Se>1.0) Se = 1.0;
	
	epsPC = 0.0001;
	
	if (Se > epsPC) 
	  return (-p0/lambda*pow(Se, -1.0/lambda-1)/(1-Snr-Swr));
	else {
	double dSedSwsquare=1/(1-Snr-Swr)/(1-Snr-Swr);
	return (-p0*dSedSwsquare/lambda/pow(epsPC, (1+1/lambda)));
	}
    }
	  
    double saturationW (double pC, double T=283.15) 
    {
      return ((1 - this->wettingPhase.Sr())*pow(p0/pC, lambda) + this->wettingPhase.Sr());
    }

    double dSdP (double pC, double T=283.15) const 
    {
      return (-(1 - this->wettingPhase.Sr())*pow(p0/pC, lambda-1)*lambda*pow(1.0/pC, 2));
    }

	  	  
    double krw (double saturationW) const
    {
      double Se = (saturationW - this->wettingPhase.Sr())
	/(1. - this->wettingPhase.Sr() - this->nonwettingPhase.Sr());
      double exponent = (2. + 3*lambda) / lambda;
      return pow(Se, exponent);
    }
      
    double krw (double saturationW, const FieldVector<double, 4>& parameters) 
    {
      double Swr = parameters[0];
      double Snr = parameters[1];
      lambda = parameters[2];
		            
      //regularisation
      //if (saturationW < (Swr + epsKr)) return 0.0;
      //if (saturationW > (1 - Snr - epsKr)) return 1.0;


      double Se = (saturationW - Swr)/(1-Snr-Swr);

      //regularisation
      if (Se > 1) return 1.0;
      if (Se < epsKr) return 0.0;

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
    double krn (double saturationN, const FieldVector<double, 4>& parameters)  
    {
      double Swr = parameters[0];
      double Snr = parameters[1];
      lambda = parameters[2];
		   
      //regularisation
      //if (saturationN < (Snr + epsKr)) return 0.0;
      //if (saturationN > (1 - Swr - epsKr)) return 1.0;

      double Se = (1-saturationN - Swr)/(1-Snr-Swr);
         
      //regularisation
      if (Se > 1) return 0.0;
      if (Se < epsKr) return 1.0;

      double exponent = (2. + lambda) / lambda;
      return pow(1-Se,2)*(1-pow(Se,exponent));
    }
		
    BrooksCoreyLaw(const Medium& wP = *(new Uniform), const Medium& nwP = *(new Uniform), 
		   double l = 2.0, double p = 0, double eps1 = 0.001, double
		   eps2=1e-15) 
      : TwoPhaseRelations(wP, nwP,
			  false),lambda(l),p0(p),epsPC(eps1),epsKr(eps2)
    {	  }

  private:
    double lambda; //!< material parameter for the relative permeability/saturation relation
    double p0; //!< capillary entry pressure
    double epsPC; //!< treshold for linearization of capillary pressure
    double epsKr;
  };
}

#endif


