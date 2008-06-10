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
      double Swr=this->wettingPhase.Sr();
      double Snr=this->nonwettingPhase.Sr();
      double pc;
      
      
      if (saturationW>(1-Snr)) return(minpc);
      if (saturationW<Swr) return(maxpc);
	    
      //pc=maxpc*(1-Snr-saturationW)/(1-Swr-Snr);
      pc=maxpc*(1-(saturationW-Swr)/(1-Swr-Snr));
 
      return pc;
    }
	
    double pC (double saturationW, const FieldVector<double, 4>& parameters) 
    {
      double Swr=this->wettingPhase.Sr();
      double Snr=this->nonwettingPhase.Sr();
      double pc;

     	    
      if (saturationW>(1-Snr)) return(minpc);
      if (saturationW<Swr) return(maxpc);
	    
      pc=maxpc*(1-(saturationW-Swr)/(1-Swr-Snr));

      return pc;
    }
	  
    double dPdS (double saturationW, double T=283.15, double p=1e5) 
    {
      double Swr=this->wettingPhase.Sr();
      double Snr=this->nonwettingPhase.Sr();
      double dpcdSw;

      dpcdSw=maxpc*(-1)/(1-Swr-Snr);

      return dpcdSw;  
    }
	  
    double saturationW (double pC, double T=283.15) 
    {
      return (0.0);
    }

    double dSdP (double pC, double T=283.15) const 
    {
      return (0.0);
    }

    LinearLaw(const Medium& wP = *(new Uniform), const Medium& nwP = *(new Uniform),double pcmax = 0.0,double pcmin = 0.0)
      : TwoPhaseRelations(wP, nwP, true),maxpc(pcmax),minpc(pcmin)
    {	 }
	  
  private:
    double maxpc;
    double minpc;
  protected:
    double krw (double saturationW) const
    { 
      return std::max(std::min((saturationW - this->wettingPhase.Sr())/(1- this->wettingPhase.Sr()-this->nonwettingPhase.Sr()), 1.), 0.);
    }
	  
    double krw (double saturationW, const FieldVector<double, 4>& parameters) 
    {
      return krw(saturationW);
    }

    double krn (double saturationN) const 
    {
      return std::max(std::min((saturationN - this->nonwettingPhase.Sr())/(1- this->wettingPhase.Sr()-this->nonwettingPhase.Sr()), 1.), 0.);
    }

    double krn (double saturationN, const FieldVector<double, 4>& parameters) 
    {
      return krn(saturationN);
    }
  };
}

#endif


