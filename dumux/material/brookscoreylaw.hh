#ifndef DUNE_BROOKSCOREYLAW_HH
#define DUNE_BROOKSCOREYLAW_HH

#include "dumux/material/relperm_pc_law.hh"

namespace Dune
{
  /*!\ingroup 2pRel 
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
   * 
   *  Vector entries in Matrix2p::paramRelPerm must be in the order
   * 		- \f$ lambda \f$
   * 		- entry pressure \f$ p_c \f$
   */
  class BrooksCoreyLaw : public RelPerm_pc
  {
  public:
    
  	double pC (double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const
    {
      double Se = (saturationW - this->soil.Sr_w(x, e, xi, T))
      	/(1. - this->soil.Sr_w(x, e, xi, T) - this->soil.Sr_n(x, e, xi, T));
      
      double lambda = this->soil.paramRelPerm(x, e, xi, T)[0];
      double p0 = this->soil.paramRelPerm(x, e, xi, T)[1];

      if (Se > epsPC) 
      	return (p0*pow(Se, -1.0/lambda));
      else 
      {
      	double dpCEps = dPdS(epsPC,283.15, 1e5);
      	return (dpCEps*(Se - epsPC) + p0*pow(epsPC, -1.0/lambda));
      }
    }
	  
    double dPdS (double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const
    {
      double Swr = this->soil.Sr_w(x, e, xi, T);
      double Snr = this->soil.Sr_n(x, e, xi, T);
      
      double lambda = this->soil.paramRelPerm(x, e, xi, T)[0];
      double p0 = this->soil.paramRelPerm(x, e, xi, T)[1];
      
      double Se = (saturationW - Swr)
      	/(1. - Swr - Snr);
      
      if(Se<0.0) Se = 0.0;
      if(Se>1.0) Se = 1.0;
	
      if (Se > epsPC) 
      	return (-p0/lambda*pow(Se, -1.0/lambda-1)/(1-Snr-Swr));
      else 
      {
      	double dSedSwsquare=1/(1-Snr-Swr)/(1-Snr-Swr);
      	return (-p0*dSedSwsquare/lambda/pow(epsPC, (1+1/lambda)));
      }
    }
	  
    double saturationW (double pC, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const
    {
    	double lambda = this->soil.paramRelPerm(x, e, xi, T)[0];
      double p0 = this->soil.paramRelPerm(x, e, xi, T)[1];
      return ((1 - this->soil.Sr_w(x, e, xi, T)) * pow(p0/pC, lambda) + this->soil.Sr_n(x, e, xi, T));
    }

    double dSdP (double pC, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const 
    {
    	double lambda = this->soil.paramRelPerm(x, e, xi, T)[0];
      double p0 = this->soil.paramRelPerm(x, e, xi, T)[1];
      return (-(1 - this->soil.Sr_w(x, e, xi, T))*pow(p0/pC, lambda-1)*lambda*pow(1.0/pC, 2));
    }

	  	  
    double krw (double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
      double Se = (saturationW - this->soil.Sr_w(x, e, xi, T))
      	/(1. - this->soil.Sr_w(x, e, xi, T) - this->soil.Sr_n(x, e, xi, T));
      
      //regularisation
      if (Se > 1) return 1.0;
      if (Se < epsKr) return 0.0;
      
      double lambda = this->soil.paramRelPerm(x, e, xi, T)[0];
      
      double exponent = (2. + 3*lambda) / lambda;
      return pow(Se, exponent);
    }
      
	  
    double krn (double saturationN, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const 
    {
      double Se = ((1.-saturationN) - this->soil.Sr_w(x, e, xi, T))
      	/(1. - this->soil.Sr_w(x, e, xi, T) - this->soil.Sr_n(x, e, xi, T));
      
      //regularisation
      if (Se > 1) return 0.0;
      if (Se < epsKr) return 1.0;
      
      double lambda = this->soil.paramRelPerm(x, e, xi, T)[0];
      double exponent = (2. + lambda) / lambda;
      return pow(1.-Se, 2) * ( 1. - pow(Se, exponent) );
    }
    
    std::vector<double> kr (const double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const
    {
    	std::vector kr(2);
    	double Srw = this->soil.Sr_w(x, e, xi, T);
    	double Srn = this->soil.Sr_n(x, e, xi, T);
    	double Se = (saturationW - Srw) / (1 - Srw - Srn);
      // regularization
    	if (Se > 1) 
  		{
    		kr[0] = 1;
    		kr[1] = 0;
    		return kr;
  		}
      if (Se < epsKr) 
    	{
      	kr[0] = 0;
    		kr[1] = 1;
    		return kr;
    	}
      
      double lambda = this->soil.paramRelPerm(x, e, xi, T)[0];
      double exponent = (2. + 3*lambda) / lambda;
      kr[0] = pow(Se, exponent);
      exponent = (2. + lambda) / lambda;
    	kr[1] = pow(1.-Se, 2) * ( 1. - pow(Se, exponent) );
    	return kr;
    }
		
    BrooksCoreyLaw(const Matrix2P& s, bool lin = false) 
      : RelPerm_pc(s, false)
    {
    }

  private:
    static const double epsPC = 0.0001;
    static const double epsKr = 1e-15;
  };
}

#endif


