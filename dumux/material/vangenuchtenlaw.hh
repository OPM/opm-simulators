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
		  double Se = std::min(1.0, (saturationW - this->wettingPhase.Sr())
				       /(1. - this->wettingPhase.Sr() - this->nonwettingPhase.Sr()));

		  if (Se > epsPC) 
			  return (p0*pow(pow(Se, -1.0/m) - 1.0, 1.0 - m));
		  else {
			  double dpCEps = dPdS(epsPC);
			  return (dpCEps*(Se - epsPC) + p0*pow(pow(epsPC, -1.0/m) - 1.0, 1.0 - m));
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
		  double Se = std::min(1.0, (saturationW - this->wettingPhase.Sr())
				       /(1. - this->wettingPhase.Sr() - this->nonwettingPhase.Sr()));
		  return (p0*(m - 1.0)/m*pow(pow(std::max(Se, epsPC), -1.0/m) - 1.0, -m)
				  *pow(std::max(Se, epsPC), -(m + 1.0)/m));
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
			  double e = 0.0, double g = 0.0, double n = 0.0, double p = 0.0, double eps = 0.0) 
			  : TwoPhaseRelations(wP, nwP, false)
	  {
	    epsilon = e ? e : 0.5;
	    gamma = g ? g : 0.5;
	    m = n ? n : 0.8;
	    p0 = p ? p : 1.0e4;
	    epsPC = eps ? eps : 5e-4;
	    machineEps_ = 1e-15;
	  }
	  
	protected:
		double krw (double saturationW) const
		{
		  double Se = std::min(1.0, (saturationW - this->wettingPhase.Sr())/(1. - this->wettingPhase.Sr() - this->nonwettingPhase.Sr()));
		    return pow( Se, epsilon ) * pow(1. - pow((1. - pow(Se, 1./m)) ,m), 2.);
		}
	  
		 double krw (double saturationW, const FieldVector<double, 4>& parameters) const 
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
		  double Se = std::min(1.0, ((1.-saturationN) - this->wettingPhase.Sr())/(1. - this->wettingPhase.Sr() - this->nonwettingPhase.Sr()));
		    return pow(1.-Se, gamma) * pow(1. - pow(Se, 1./m), 2.*m);
		}
		
		 double krn (double saturationN, const FieldVector<double, 4>& parameters) const 
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
		double m; //!< parameter for the relative permeability/saturation relation
		double p0; //!< capillary entry pressure
		double epsPC; //!< treshold for linearization of capillary pressure
		double machineEps_;
	};
}

#endif


