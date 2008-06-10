#ifndef DUNE_VANGENUCHTENLAW_HH
#define DUNE_VANGENUCHTENLAW_HH

#include "dumux/material/twophaserelations.hh"

namespace Dune
{
  /*!\ingroup material 
   * \brief van Genuchten mobility/saturation relation.  
   *
   *  Employs the van Genuchten non-linear relative permeability/saturation relation.
   */
  class VanGenuchtenLaw : public TwoPhaseRelations
  {
  public:
    double pC (double saturationW) 
    {
      double r,x,vgM;
      double pc,pc_prime,Se_regu;
      int asymptotic;

      double Swr = this->wettingPhase.Sr();
      double Snr = this->nonwettingPhase.Sr();

      double Se = (saturationW - Swr)/(1. - Swr - Snr);

      if(Se<0.0) Se = 0.0;
      if(Se>1.0) Se = 1.0;

      /* check size of S_e^(-1/m) */
      if (Se<epsPC)
	asymptotic=1;
      else if (Se>1-epsPC)
	asymptotic=0;
      else if (1.0/pow(Se,1.0/m)>1000.0)
	asymptotic=1;
      else asymptotic=0;

      if (asymptotic) /* use pc_VanG = 1/alpha/pow(Se,1/(n-1)) */
	{
	  if (Se>epsPC)
	    return(1.0/(alpha*pow(Se,1.0/(n-1))));
	  else /* regularize with tangent */
	    {
	      pc  = 1.0/(alpha*pow(epsPC,1.0/(n-1)));
	      pc_prime = 1.0/(pow(epsPC,n/(n-1))*alpha*(1-n)*(1-Swr-Snr));
	      return((Se-epsPC)*pc_prime+pc);
	    }
	}
      else
	{	/* use correct van Genuchten curve */
	  if (Se>epsPC && Se<1-epsPC)
	    {
	      r = pow(Se,-1/m);
	      x = r-1;
	      vgM = 1-m;
	      x = pow(x,vgM);
	      r = x/alpha;
	      return(r);
	    }
	  else
	    {
	      /* value and derivative at regularization point */
	      if (Se<=epsPC) Se_regu = epsPC; else Se_regu = 1-epsPC;
	      pc       = pow(pow(Se_regu,-1/m)-1,1/n)/alpha;
	      pc_prime = pow(pow(Se_regu,-1/m)-1,1/n-1)*pow(Se_regu,-1/m-1)*(-1/m)/alpha/(1-Snr-Swr);
			
	      /* evaluate tangential */
	      r        = (Se-Se_regu)*pc_prime+pc;
	      return(r);
	    }
	}
    }
	
    double pC (double saturationW, const FieldVector<double, 4>& parameters) 
    {
      double Swr = parameters[0];
      double Snr = parameters[1];
      double alpha = parameters[2];
      double vgN = parameters[3];
		  
      double r,Se,x,vgM;
      double pc,pc_prime,Se_regu;
      int asymptotic;
			
      /* effective values */
      Se = (saturationW-Swr)/(1-Snr-Swr);
      if(Se<0.0) Se = 0.0;
      if(Se>1.0) Se = 1.0;
      vgM=1-1/vgN;
			
      /* check size of S_e^(-1/m) */
      if (Se<epsPC)
	asymptotic=1;
      else if (Se>1-epsPC)
	asymptotic=0;
      else if (1.0/pow(Se,1.0/vgM)>1000.0)
	asymptotic=1;
      else asymptotic=0;

      if (asymptotic) /* use pc_VanG = 1/alpha/pow(Se,1/(n-1)) */
	{
	  if (Se>epsPC)
	    return(1.0/(alpha*pow(Se,1.0/(vgN-1))));
	  else /* regularize with tangent */
	    {
	      pc       = 1.0/(alpha*pow(epsPC,1.0/(vgN-1)));
	      pc_prime = 1.0/(pow(epsPC,vgN/(vgN-1))*alpha*(1-vgN)*(1-Swr-Snr));
	      return((Se-epsPC)*pc_prime+pc);
	    }
	}
      else
	{	/* use correct van Genuchten curve */
	  if (Se>epsPC && Se<1-epsPC)
	    {
	      r = pow(Se,-1/vgM);
	      x = r-1;
	      vgM = 1-vgM;
	      x = pow(x,vgM);
	      r = x/alpha;
	      return(r);
	    }
	  else
	    {
	      /* value and derivative at regularization point */
	      if (Se<=epsPC) Se_regu = epsPC; else Se_regu = 1-epsPC;
	      pc       = pow(pow(Se_regu,-1/vgM)-1,1/vgN)/alpha;
	      pc_prime = pow(pow(Se_regu,-1/vgM)-1,1/vgN-1)*pow(Se_regu,-1/vgM-1)*(-1/vgM)/alpha/(1-Snr-Swr);
			
	      /* evaluate tangential */
	      r        = (Se-Se_regu)*pc_prime+pc;
	      return(r);
	    }
	}
    }
	
    double dPdS (double saturationW, double T=283.15, double p=1e5)
    {
      double r,x;

      double Swr = this->wettingPhase.Sr();
      double Snr = this->nonwettingPhase.Sr();
      
      double Se = (saturationW - Swr)/(1. -Swr - Snr);
      
      /* regularization */
      if (Se<=0.0) Se=1.E-3;
      if (Se>=1.0) Se=1.0-1.E-5;

      /* compute value */
        r = pow(Se,-1/m);
        x = pow((r-1),(1-m));
        r = -(1-0.0)/alpha*x*(1-m)/m*r/(r-1)/Se/(1-Snr-Swr);
        return(r);
    }
	  
    double saturationW (double pC, double T=283.15) 
    {
      DUNE_THROW(NotImplemented, "VanGenuchtenLaw::saturationW()");
		  
      return 0;
    }

    double dSdP (double pC, double T=283.15) const 
    {
      DUNE_THROW(NotImplemented, "VanGenuchtenLaw::dSdP()");
		  
      return 0;
    }

    VanGenuchtenLaw(const Medium& wP = *(new Uniform), const Medium& nwP = *(new Uniform), 
		    double n = 0.0, double a = 0.0, double e = 0.0, double g = 0.0, double eps = 0.0) 
      : TwoPhaseRelations(wP, nwP, false)
    {
      epsilon = e ? e : 0.5;
      gamma = g ? g : (1/3);
      n = n ? n : 4.4;
      m = (1+1/n);
      alpha = a ? a : 0.4;
      epsPC = eps ? eps : 5e-4;
      machineEps_ = 1e-15;
    }
	  
  protected:
    double krw (double saturationW) const
    {
      double Swr = this->wettingPhase.Sr();
      double Snr = this->nonwettingPhase.Sr();

      double Se,krw,r;

      /* regularization */
      if(saturationW<(Swr+machineEps_)) r=0.;
      if(saturationW>(1.-Snr-machineEps_)) r=1.;

      /* effective values */
      Se = (saturationW-Swr)/(1-Snr-Swr);               

      /* regularization */
      if(Se>1.)  Se=1.;
      if(Se<machineEps_) Se=machineEps_;

      /* compute value */
      r   = 1-pow(Se,1/m);
      krw = sqrt(Se)*pow(1-pow(r,m),2);
      return(krw);
    }
	  
    double krw (double saturationW, const FieldVector<double, 4>& parameters)  
    {
      double Swr = parameters[0];
      double Snr = parameters[1];
      //double alpha = parameters[2];
      double vgN = parameters[3];

      double Se,vgM,krw,r;

      /* regularization */
      if(saturationW<(Swr+machineEps_)) r=0.;
      if(saturationW>(1.-Snr-machineEps_)) r=1.;

      /* effective values */
      Se = (saturationW-Swr)/(1-Snr-Swr);
      vgM = 1-1/vgN;                 /* >3 */

      /* regularization */
      if(Se>1.)  Se=1.;
      if(Se<machineEps_) Se=machineEps_;

      /* compute value */
      r   = 1-pow(Se,1/vgM);
      krw = sqrt(Se)*pow(1-pow(r,vgM),2);
      return(krw);
    }

    double krn (double saturationN) const 
    {
      double Swr = this->wettingPhase.Sr();
      double Snr = this->nonwettingPhase.Sr();

      double Se,krn,r;

      /* regularization */
      if(saturationN<(Snr+machineEps_)) r=0.;
      if(saturationN>(1.-Swr-machineEps_)) r=1.;

      /* effective values */
      Se = (1-saturationN-Swr)/(1-Snr-Swr);

      /* effective Saturation Se has to be between 0 and 1! */
      if(Se>1.)  Se=1.;
      if(Se<machineEps_) Se=machineEps_;

      /* compute value */
      r   = 1-pow(Se,1/m);
      krn = sqrt(1-Se)*pow(pow(r,m),2);
      return(krn);
    }
		
    double krn (double saturationN, const FieldVector<double, 4>& parameters)  
    {
      double Swr = parameters[0];
      double Snr = parameters[1];
      //double alpha = parameters[2];
      double vgN = parameters[3];
			  
      double Se,vgM,krn,r;

      /* regularization */
      if(saturationN<(Snr+machineEps_)) r=0.;
      if(saturationN>(1.-Swr-machineEps_)) r=1.;

      /* effective values */
      Se = (1-saturationN-Swr)/(1-Snr-Swr);
      vgM = 1-1/vgN;                /* >3 */

      /* effective Saturation Se has to be between 0 and 1! */
      if(Se>1.)  Se=1.;
      if(Se<machineEps_) Se=machineEps_;

      /* compute value */
      r   = 1-pow(Se,1/vgM);
      krn = sqrt(1-Se)*pow(pow(r,vgM),2);
      return(krn);
    }
  
    double epsilon; //!< material parameter for the wetting phase
    double gamma; //!< material parameter for the non-wetting phase
    double n;
    double m; //!< parameter for the relative permeability/saturation relation
    double alpha; //!< capillary entry pressure
    double epsPC; //!< treshold for linearization of capillary pressure
    double machineEps_;
  };
}

#endif


