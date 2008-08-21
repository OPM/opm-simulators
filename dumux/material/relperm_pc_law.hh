#ifndef DUNE_RELATIVEPERMEABILITY_HH
#define DUNE_RELATIVEPERMEABILITY_HH

#include <dumux/material/property_baseclasses.hh>

/**
 * \author Jochen Fritz
 */


namespace Dune
{
	/*!\ingroup 2pRel 
	 * \brief Base class for relative permeability and capillary pressure - saturation relationships.  
	 * Derived classes of this base class are used as input for the TwoPhaseRelations class.
	 * The functions get significant input from the member soil which is an object of Matrix2p.
	 * Specifications of model parameters for the relative permeability and capillary pressure
	 * functions have to be made in right order. Read further details in the derived classes!
	 */
	template<class G>
	class RelPerm_pc {
	public:
		typedef typename G::Traits::template Codim<0>::Entity Entity;
		typedef typename G::ctype DT;
		enum {n=G::dimension, m=1};
	  /*! \brief the capillary pressure - saturation relation 
	   *
	   *  \param saturationW the saturation of the wetting phase  
	   *  \param x position in global coordinates
	   *  \param e codim 0 entity for which the value is sought
	   *  \param xi position in local coordinates in e 
	   *  \return the capillary pressur \f$ p_\text{c} (S_\text{w})\f$.
	   */
	  virtual double pC (double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T) const = 0;

	  /*! \brief the derivative of capillary pressure w.r.t. the saturation 
	   *
	   *  \param saturationW the saturation of the wetting phase   
	   *  \param T temperature
	   *  \param p pressure   
	   *  \param x position in global coordinates
	   *  \param e codim 0 entity for which the value is sought
	   *  \param xi position in local coordinates in e 
	   *  \return the derivative \f$\text{d}p_\text{c}/\text{d}S_\text{e}\f$
	   */
	  virtual double dPdS (double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T) const = 0;
	  
	  /*! \brief the wetting phase saturation w.r.t. the capillary pressure 
	   *
	   *  \param pC the capillary pressure   
	   *  \param T temperature
	   *  \param x position in global coordinates
	   *  \param e codim 0 entity for which the value is sought
	   *  \param xi position in local coordinates in e 
	   *  \return the wetting phase saturation 
	   */
	  virtual double saturationW (double pC, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T) const = 0;
	  
	  /*! \brief the derivative of the saturation w.r.t. the capillary pressure 
	   *
	   *  \param pC the capillary pressure   
	   *  \param T temperature
	   *  \param x position in global coordinates
	   *  \param e codim 0 entity for which the value is sought
	   *  \param xi position in local coordinates in e 
	   *  \return the derivative \f$\text{d}S_w/\text{d}p_\text{c}\f$
	   */
	  virtual double dSdP (double pC, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T) const = 0;
	  
	  const bool isLinear() const 
	  {
		  return linear_;
	  }
	  
		/*! \brief wetting phase relative permeability saturation relationship 
		 *
		 *  \param saturationW the saturation of the wetting phase
	   *  \param x position in global coordinates
	   *  \param e codim 0 entity for which the value is sought
	   *  \param xi position in local coordinates in e    
		 *  \return the wetting phase relative permeability
		 */
		virtual double krw (const double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const = 0; 
	  
		/*! \brief nonwetting phase relative permeability saturation relationship 
		 *
		 *  \param saturationN the saturation of the nonwetting phase
	   *  \param x position in global coordinates
	   *  \param e codim 0 entity for which the value is sought
	   *  \param xi position in local coordinates in e    
		 *  \return the nonwetting phase relative permeability
		 */
		virtual double krn (const double saturationN, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const = 0;
	  
		
		/** \brief relative permeability saturation relationship for both phases
		 *  In many cases the relative permeabilities of both phases are needed at 
		 *  the same time. This function reduces unnecessary computational costs.
		 *  \param saturationN the saturation of the nonwetting phase
	   *  \param x position in global coordinates
	   *  \param e codim 0 entity for which the value is sought
	   *  \param xi position in local coordinates in e    
		 *  \return relative permeability vector: first entry wetting, seconde entry nonwetting phase
		 */
		virtual std::vector<double> kr (const double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const = 0;
		
		/** \brief constructor
		 *  \param s a matrix property object.
		 *  \param wP phase property object for the wetting Phase.
		 *  \param nP phase property object for the nonwetting Phase.
		 *  \param lin true specifies a linear model. Usually false. Only set true if you know what you are doing!
		 */
	  RelPerm_pc(const Matrix2p<G,double>& s, const bool lin = false)
	  : soil(s), linear_(lin)
	  {	 
	  }
	  
	  virtual ~RelPerm_pc()
	  {	  
	  }
	
	protected:
		const bool linear_;
	  const Matrix2p<G,double>& soil;
	};
}

#endif
