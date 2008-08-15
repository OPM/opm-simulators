#ifndef DUNE_LINEARLAW_HH
#define DUNE_LINEARLAW_HH

#include "dumux/material/twophaserelations.hh"

namespace Dune
{
  /** \ingroup 2pRel
   * @brief Represents the linear relative permability - saturation relation \f$ k_{r\alpha}(S) = \frac{S-S_{r\alpha}}{1-S_{r\alpha}} \f$
   * Vector entries in Matrix2p::paramRelPerm must be in the order
   * 		- minimum capillary pressure
   * 		- maximum capillary pressure  
   */
  class LinearLaw : public TwoPhaseRelations 
  {
  public:
  	
    double krw (double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const
    { 
    	double Sr_w = this->soil.Sr_w(x, e, xi, T);
    	double Sr_n = this->soil.Sr_n(x, e, xi, T);
      return std::max(std::min((saturationW - Sr_w)/(1- Sr_w - Sr_n), 1.), 0.);
    }

    double krn (double saturationN, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const 
    {
    	double Sr_w = this->soil.Sr_w(x, e, xi, T);
    	double Sr_n = this->soil.Sr_n(x, e, xi, T);
      return std::max(std::min((saturationN - Sr_n)/(1- Sr_w - Sr_n), 1.), 0.);
    }
  	
    std::vector<double> kr (const double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const
    {
    	std::vector kr(2);
    	double Sr_w = this->soil.Sr_w(x, e, xi, T);
    	double Sr_n = this->soil.Sr_n(x, e, xi, T);
    	kr[0] = std::max(std::min((saturationW - Sr_w)/(1- Sr_w - Sr_n), 1.), 0.);
    	kr[1] = std::max(std::min((saturationN - Sr_n)/(1- Sr_w - Sr_n), 1.), 0.);
    	return kr;
    }
    
    double pC (double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) 
    {
      double Swr = this->soil.Sr_w(x, e, xi, T);
      double Snr = this->soil.Sr_n(x, e, xi, T);      
      
      std::vector<double> param = this->soil.paramRelPerm(x, e, xi, T);
      
      if (saturationW > (1-Snr)) return param[0]; // min pc
      if (saturationW < Swr) return param[1]; // max pc
	    
      return pc = param[0] + (param[1] - param[0]) * (1 - (saturationW-Swr) / (1-Swr-Snr));
    }
	  
    double dPdS (double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) 
    {
      double Swr = this->soil.Sr_w(x, e, xi, T);
      double Snr = this->soil.Sr_n(x, e, xi, T); 

      std::vector<double> param = this->soil.paramRelPerm(x, e, xi, T);
      
      return dpcdSw = (param[1] - param[0]) * (-1)/(1-Swr-Snr);
    }
	  
    double saturationW (double pC, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) 
    {
    	DUNE_THROW(NotImplemented, "LinearLaw::saturationW");
    }

    double dSdP (double pC, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const 
    {
    	DUNE_THROW(NotImplemented, "LinearLaw::dSdP");
    }

    LinearLaw(const Matrix2P& s, bool lin = false) 
    : TwoPhaseRelations(s, false)
    {	 }
	  
  private:
    double maxpc;
    double minpc;
  };
}

#endif


