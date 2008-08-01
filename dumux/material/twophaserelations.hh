#ifndef DUNE_RELATIVEPERMEABILITY_HH
#define DUNE_RELATIVEPERMEABILITY_HH

#include<dumux/material/properties.hh>

/**
 * \defgroup material Relative Permeabilities and material properties
 */


namespace Dune
{
	/*!\ingroup material 
	 * \brief Base class for relative permeability - saturation relationships.  
	 *
	 *  Depending on the saturation \f$S_\text{w}\f$ and the viscosities \f$\mu_\text{w}\f$,\f$\mu_\text{n}\f$, 
	 *  the phase mobilities \f$\lambda_\text{w}\f$, \f$\lambda_\text{n}\f$, 
	 *  the total mobility \f$\lambda\f$, the capillary pressure \f$p_\text{c}\f$, and its derivative 
	 *  \f$\text{d}p_\text{c}/\text{d}S_\text{e}\f$ can be calculated. 
	 *  For use in Properties.
	 */
	class TwoPhaseRelations {
	public:
	  /*! \brief Implements the wetting phase mobility/saturation relation. 
	   *
	   *  \param saturationW the saturation of the wetting phase   
	   *  \param T temperature
	   *  \param p pressure   
	   *  \return the mobility of the wetting phase
	   */
	  double mobW (double saturationW, double T=283.15, double p=1e5)
	  {
	    double viscosityW = wettingPhase.viscosity ( T, p );
	    return krw(saturationW)/viscosityW;
	  } 

	  double mobW (double saturationW, const FieldVector<double, 4>& parameters, double T=283.15, double p=1e5)
	  {
	    double viscosityW = wettingPhase.viscosity ( T, p );
	    return krw(saturationW, parameters)/viscosityW;
	  } 

	  /*! \brief Implements the nonwetting phase mobility/saturation relation. 
	   *
	   *  \param saturationN the saturation of the nonwetting phase   
	   *  \param T temperature
	   *  \param p pressure   
	   *  \return the mobility of the nonwetting phase
	   */
	  double mobN (double saturationN, double T=283.15, double p=1e5)
	  {
	    double viscosityN = nonwettingPhase.viscosity ( T, p );
	    return krn(saturationN)/viscosityN;
	  }
	  
	  double mobCO2 (double saturationN, double T=283.15, double p=1e5, double rho=500.)
	  {
	    double viscosityN = nonwettingPhase.viscosity ( T, p, rho);
	    return krn(saturationN)/viscosityN;
	  }
	  
	  double mobN (double saturationN, const FieldVector<double, 4>& parameters, double T=283.15, double p=1e5)
	  {
	    double viscosityN = nonwettingPhase.viscosity ( T, p );
	    return krn(saturationN, parameters)/viscosityN;
	  } 
	  
	  double mobCO2 (double saturationN, const FieldVector<double, 4>& parameters, double T=283.15, double p=1e5, double rho=500.)
	  {
	    double viscosityN = nonwettingPhase.viscosity ( T, p, rho );
	    return krn(saturationN, parameters)/viscosityN;
	  } 
	  /*! \brief Implements the total mobility/saturation relation. 
	   *
	   *  \param saturationW the saturation of the wetting phase   
	   *  \param T temperature
	   *  \param p pressure   
	   *  \return the total mobility 
	   */
	  double mobTotal (double saturationW, double T=283.15, double p=1e5) 
	  {
		  return mobW(saturationW, T, p) + mobN(1.-saturationW, T, p);    
	  }
	
     double mobTotal (double saturationW, const FieldVector<double, 4>& parameters, double T=283.15, double p=1e5) 
	  {
		  return mobW(saturationW, parameters, T, p) + mobN(1.-saturationW, parameters, T, p);    
	  }

	  /*! \brief Implements the wetting phase fractional flow function
	   *
	   *  \param saturationW the saturation of the wetting phase   
	   *  \param T temperature
	   *  \param p pressure   
	   *  \return the fractional flow of the wetting phase
	   */
	  double fractionalW (double saturationW, double T=283.15, double p=1e5)
	  { 
	    double l_w = mobW(saturationW, T, p); 
	    double l_n = mobN(1.-saturationW, T, p);
	    return l_w/(l_w + l_n);
	  }
	 
     double fractionalW (double saturationW, const FieldVector<double, 4>& parameters, double T=283.15, double p=1e5)
	  { 
	    double l_w = mobW(saturationW, parameters, T, p); 
	    double l_n = mobN(1.-saturationW, parameters, T, p);
	    return l_w/(l_w + l_n);
	  }

	  /*! \brief Implements the nonwetting phase fractional flow function
	   *
	   *  \param saturationN the saturation of the nonwetting phase   
	   *  \param T temperature
	   *  \param p pressure   
	   *  \return the fractional flow of the nonwetting phase
	   */
	  double fractionalN (double saturationN, double T=283.15, double p=1e5)
	  { 
	    double l_w = mobW(1.-saturationN, T, p); 
	    double l_n = mobN(saturationN, T, p);
	    return l_n/(l_w + l_n);
	  }
	
     double fractionalN (double saturationN, const FieldVector<double, 4>& parameters, double T=283.15, double p=1e5)
	  { 
	    double l_w = mobW(1.-saturationN, parameters, T, p); 
	    double l_n = mobN(saturationN, parameters, T, p);
	    return l_n/(l_w + l_n);
	  }

	  /*! \brief the capillary pressure - saturation relation 
	   *
	   *  \param saturationW the saturation of the wetting phase   
	   *  \return the capillary pressur \f$ p_\text{c} (S_\text{w})\f$.
	   */
	  virtual double pC (double saturationW) = 0;
	
	  virtual double pC (double saturationW, const FieldVector<double, 4>& parameters) 
	  {
		  DUNE_THROW(NotImplemented,"TwoPhaseRelations :: pC (double, const FieldVector<double, 4>&");
	  }

	  /*! \brief the derivative of capillary pressure w.r.t. the saturation 
	   *
	   *  \param saturationW the saturation of the wetting phase   
	   *  \param T temperature
	   *  \param p pressure   
	   *  \return the derivative \f$\text{d}p_\text{c}/\text{d}S_\text{e}\f$
	   */
	  virtual double dPdS (double saturationW, double T=283.15, double p=1e5) = 0;
	  
	  /*! \brief the wetting phase saturation w.r.t. the capillary pressure 
	   *
	   *  \param pC the capillary pressure   
	   *  \param T temperature
	   *  \return the wetting phase saturation 
	   */
	  virtual double saturationW (double pC, double T=283.15) = 0;
	  
	  /*! \brief the derivative of the saturation w.r.t. the capillary pressure 
	   *
	   *  \param pC the capillary pressure   
	   *  \param T temperature
	   *  \return the derivative \f$\text{d}S_w/\text{d}p_\text{c}\f$
	   */
	  virtual double dSdP (double pC, double T=283.15) const = 0;
	  
	  const bool isLinear() const 
	  {
		  return linear_;
	  }
	  
	  TwoPhaseRelations(const Medium& wP = *(new Uniform), const Medium& nwP = *(new Uniform), 
			  			const bool lin = false)
	  : wettingPhase(wP), nonwettingPhase(nwP), linear_(lin)
	  {	 }
	  
	  virtual ~TwoPhaseRelations()
	  {	  }
	  
		
	  const Medium& wettingPhase; //!< contains the properties of the wetting phase 
	  const Medium& nonwettingPhase; //!< contains the properties of the nonwetting phase 
	protected:
		double p0; //!< capillary entry pressure
		const bool linear_;
		
		/*! \brief wetting phase relative permeability saturation relationship 
		 *
		 *  \param saturationW the saturation of the wetting phase   
		 *  \return the wetting phase relative permeability
		 */
		virtual double krw (double saturationW) const=0; 
	  
		virtual double krw (double saturationW, const FieldVector<double, 4>& parameters)
		  {
			  DUNE_THROW(NotImplemented,"TwoPhaseRelations :: krw (double, const FieldVector<double, 4>&");
		  }
	  
		/*! \brief nonwetting phase relative permeability saturation relationship 
		 *
		 *  \param saturationN the saturation of the nonwetting phase   
		 *  \return the nonwetting phase relative permeability
		 */
		virtual double krn (double saturationN) const=0;

		virtual double krn (double saturationN, const FieldVector<double, 4>& parameters)
		  {
			  DUNE_THROW(NotImplemented,"TwoPhaseRelations :: krn (double, const FieldVector<double, 4>&");
		  }
	};
}

#endif
