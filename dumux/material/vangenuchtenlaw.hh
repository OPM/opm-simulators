#ifndef DUNE_VANGENUCHTENLAW_HH
#define DUNE_VANGENUCHTENLAW_HH

#include "dumux/material/twophaserelations.hh"

namespace Dune
{
  /*!\ingroup 2pRel 
   * \brief van Genuchten mobility/saturation relation.  
   *
   *  Employs the van Genuchten non-linear relative permeability/saturation relation.
   *  Vector entries in Matrix2p::paramRelPerm must be in the order
   * 		- m
   * 		- n
   * 		- \f$ \epsilon \f$
   * 		- \f$ \gamma \f$
   * 		- \f$ \alpha \f$
   *   
   */
  class VanGenuchtenLaw : public TwoPhaseRelations
  {
  public:
  	
    double krw (double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const
    {
      double Swr = this->soil.Sr_w(x, e, xi, T);
      double Snr = this->soil.Sr_n(x, e, xi, T);

      double Se,krw,r;

      /* effective values */
      Se = (saturationW - Swr) / (1 - Snr - Swr);               

      /* regularization */
      if(Se > 1.) return 1.;
      if(Se < machineEps_) Se = machineEps_;

      std::vector<double> param = this->soil.paramRelPerm(x, e, xi, T);
      double m = param[0];
      double eps = param[2];
      
      /* compute value */
      r   = 1 - pow(Se, 1/m);
      krw = pow(Se, eps) * pow(1 - pow(r,m), 2);
      return(krw);
    }

    double krn (double saturationN, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const 
    {
      double Swr = this->soil.Sr_w(x, e, xi, T);
      double Snr = this->soil.Sr_n(x, e, xi, T);

      double Se, r;

      /* effective values */
      Se = (1 - saturationN - Swr) / (1 - Snr - Swr);

      /* effective Saturation Se has to be between 0 and 1! */
      if(Se > 1.) Se = 1.;
      if(Se < machineEps_) Se = machineEps_;

      std::vector<double> param = this->soil.paramRelPerm(x, e, xi, T);
      double m = param[0];
      double gamma = param[3];
      
      /* compute value */
      r   = 1 - pow(Se, 1/m);
      return pow(1-Se, gamma) * pow(r, 2*m);
    }
    
    std::vector<double> kr (const double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const
    {
    	std::vector kr(2);
    	// residual saturations
   		double Srw = this->soil.Sr_w(x, e, xi, T);
   		double Srn = this->soil.Sr_n(x, e, xi, T);
   		// effective saturation
   		double Se = (saturationW - Srw) / (1 - Srw - Srn);
   		
   		// regularization
      if(Se > 1.)
      {
      	kr[0] = 1;
      	kr[1] = 0;
      }
      if(Se < machineEps_) Se = machineEps_;
      
      // get Van Genuchten parameters
      std::vector<double> param = this->soil.paramRelPerm(x, e, xi, T);
      double m = param[0];
      double eps = param[2];
      double gamma = param[3];
      
      // compute values
      r   = 1 - pow(Se, 1/m);
      kr[0] = pow(Se, eps) * pow(1 - pow(r,m), 2);
      kr[1] = pow(1-Se, gamma) * pow(r, 2*m);
      return kr;
    }
    
  	double pC (const double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const 
    {
      double r, x, vgM;
      double pc, pc_prime, Se_regu;
      int asymptotic;

      double Swr = this->soil.Sr_w(x, e, xi, T);
      double Snr = this->soil.Sr_n(x, e, xi, T);
      
      std::vector<double> param = this->soil.paramRelPerm(x, e, xi, T);
      double m = param[0];
      double n = param[1];
      double alpha = param[5];

      double Se = (saturationW - Swr) / (1. - Swr - Snr);

      if(Se < 0.0) Se = 0.0;
      if(Se > 1.0) Se = 1.0;

      /* check size of S_e^(-1/m) */
      if (Se < epsPC)
      	asymptotic = 1;
      else if (Se > 1 - epsPC)
      	asymptotic = 0;
      else if (1.0 / pow(Se, 1.0/m) > 1000.0)
      	asymptotic = 1;
      else asymptotic = 0;

      if (asymptotic) /* use pc_VanG = 1/alpha/pow(Se,1/(n-1)) */
      {
      	if (Se > epsPC)
      		return(1.0 / (alpha * pow(Se, 1.0 / (n-1))));
      	else /* regularize with tangent */
      	{
		      pc  = 1.0 / (alpha * pow(epsPC, 1.0 / (n-1)));
		      pc_prime = 1.0 / (pow(epsPC, n/(n-1)) * alpha * (1-n) * (1-Swr-Snr));
		      return((Se - epsPC) * pc_prime + pc);
      	}
      }
      else
      {	/* use correct van Genuchten curve */
      	if (Se > epsPC && Se < 1 - epsPC)
      	{
		      r = pow(Se, -1/m);
		      x = r - 1;
		      vgM = 1 - m;
		      x = pow(x, vgM);
		      r = x / alpha;
		      return(r);
      	}
      	else
      	{
		      /* value and derivative at regularization point */
		      if (Se <= epsPC) Se_regu = epsPC; else Se_regu = 1 - epsPC;
		      pc       = pow(pow(Se_regu, -1/m) - 1, 1/n) / alpha;
		      pc_prime = pow(pow(Se_regu, -1/m) - 1, 1/n-1) * pow(Se_regu, -1/m-1) * (-1/m) / alpha / (1-Snr-Swr);
				
		      /* evaluate tangential */
		      r        = (Se - Se_regu) * pc_prime + pc;
		      return(r);
      	}
      }
    }
	
    
  	double dPdS (double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const
    {
      double r, x;
s
      double Swr = this->soil.Sr_w(x, e, xi, T);
      double Snr = this->soil.Sr_n(x, e, xi, T);
      
      std::vector<double> param = this->soil.paramRelPerm(x, e, xi, T);
      double m = param[0];
      double n = param[1];
      double alpha = param[5];
      
      double Se = (saturationW - Swr) / (1. -Swr - Snr);
      
      /* regularization */
      if (Se <= 0.0) Se = 1.E-3;
      if (Se >= 1.0) Se = 1.0 - 1.E-5;

      /* compute value */
        r = pow(Se, -1/m);
        x = pow((r-1), (1-m));
        r = -(1-0.0) / alpha * x * (1-m) / m * r / (r-1) / Se / (1-Snr-Swr);
        return(r);
    }
	  
    double saturationW (double pC, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) 
    {
      DUNE_THROW(NotImplemented, "VanGenuchtenLaw::saturationW()");
		  
      return 0;
    }

    double dSdP (double pC, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const 
    {
      DUNE_THROW(NotImplemented, "VanGenuchtenLaw::dSdP()");
		  
      return 0;
    }

    VanGenuchtenLaw(const Matrix2P& s, bool lin = false) 
      : TwoPhaseRelations(s, false)
    {
    }
	  
  protected:
    double epsPC = 5e-4; //!< threshold for linearization of capillary pressure
    double machineEps_ = 1e-15;
  };
}

#endif


