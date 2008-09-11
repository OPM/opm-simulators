// $Id$ 

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


	 /** \ingroup 2pRel
	 * @brief Represents the linear relative permability - saturation relation \f$ k_{r\alpha}(S) = \frac{S-S_{r\alpha}}{1-S_{r\alpha}} \f$
	 * Vector entries in Matrix2p::paramRelPerm must be in the order
	 * 		- minimum capillary pressure
	 * 		- maximum capillary pressure
	 */
	template<class G>
	class LinearLaw : public RelPerm_pc<G>
	{
	public:
		typedef typename G::Traits::template Codim<0>::Entity Entity;
		typedef typename G::ctype DT;
		enum {n=G::dimension, m=1};

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
			std::vector<double> kr(2);
			double Sr_w = this->soil.Sr_w(x, e, xi, T);
			double Sr_n = this->soil.Sr_n(x, e, xi, T);
			kr[0] = std::max(std::min((saturationW - Sr_w)/(1- Sr_w - Sr_n), 1.), 0.);
			kr[1] = std::max(std::min((1 - saturationW - Sr_n)/(1- Sr_w - Sr_n), 1.), 0.);
			return kr;
		}

		double pC (double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const
		{
			double Swr = this->soil.Sr_w(x, e, xi, T);
			double Snr = this->soil.Sr_n(x, e, xi, T);

			std::vector<double> param = this->soil.paramRelPerm(x, e, xi, T);

			if (saturationW > (1-Snr)) return param[0]; // min pc
			if (saturationW < Swr) return param[1]; // max pc

			return  param[0] + (param[1] - param[0]) * (1 - (saturationW-Swr) / (1-Swr-Snr));
		}

		double dPdS (double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const
		{
			double Swr = this->soil.Sr_w(x, e, xi, T);
			double Snr = this->soil.Sr_n(x, e, xi, T);

			std::vector<double> param = this->soil.paramRelPerm(x, e, xi, T);

			return (param[1] - param[0]) * (-1)/(1-Swr-Snr);
		}

		double saturationW (double pC, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const
		{
			DUNE_THROW(NotImplemented, "LinearLaw::saturationW");
		}

		double dSdP (double pC, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const
		{
			DUNE_THROW(NotImplemented, "LinearLaw::dSdP");
		}

		LinearLaw(const Matrix2p<G,double>& s, bool lin = false)
		: RelPerm_pc<G>(s, false)
		{	 }

	private:
		double maxpc;
		double minpc;
	};


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
	template<class G>
	class VanGenuchtenLaw : public RelPerm_pc<G>
	{
	public:

		typedef typename G::Traits::template Codim<0>::Entity Entity;
		typedef double DT;
		enum {n=G::dimension, m=1};

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
			std::vector<double> kr(2);
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
			double r   = 1 - pow(Se, 1/m);
			kr[0] = pow(Se, eps) * pow(1 - pow(r,m), 2);
			kr[1] = pow(1-Se, gamma) * pow(r, 2*m);
			return kr;
		}

		double pC (const double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const
		{
			double r, x_, vgM;
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
					x_ = r - 1;
					vgM = 1 - m;
					x_ = pow(x_, vgM);
					r = x_ / alpha;
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
			double r, x_;

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
				x_ = pow((r-1), (1-m));
				r = -(1-0.0) / alpha * x_ * (1-m) / m * r / (r-1) / Se / (1-Snr-Swr);
				return(r);
		}

		double saturationW (double pC, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const
		{
			DUNE_THROW(NotImplemented, "VanGenuchtenLaw::saturationW()");

			return 0;
		}

		double dSdP (double pC, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, double T) const
		{
			DUNE_THROW(NotImplemented, "VanGenuchtenLaw::dSdP()");

			return 0;
		}

		VanGenuchtenLaw(const Matrix2p<G,double>& s, bool lin = false)
			: RelPerm_pc<G>(s, false)
		{
		}

	protected:
		static const double epsPC = 5e-4; //!< threshold for linearization of capillary pressure
		static const double machineEps_ = 1e-15;
	};


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
	template<class G>
	class BrooksCoreyLaw : public RelPerm_pc<G>
	{
	public:
		typedef typename G::Traits::template Codim<0>::Entity Entity;
		typedef typename G::ctype DT;
		enum {n=G::dimension, m=1};

		double pC (double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T) const
		{
			double Se = (saturationW - this->soil.Sr_w(x, e, xi, T))
				/(1. - this->soil.Sr_w(x, e, xi, T) - this->soil.Sr_n(x, e, xi, T));

			double lambda = this->soil.paramRelPerm(x, e, xi, T)[0];
			double p0 = this->soil.paramRelPerm(x, e, xi, T)[1];

			if (Se > epsPC)
				return (p0*pow(Se, -1.0/lambda));
			else
			{
				double dpCEps = dPdS(epsPC, x, e, xi, T);
				return (dpCEps*(Se - epsPC) + p0*pow(epsPC, -1.0/lambda));
			}
		}

		double dPdS (double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T) const
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

		double saturationW (double pC, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T) const
		{
			double lambda = this->soil.paramRelPerm(x, e, xi, T)[0];
			double p0 = this->soil.paramRelPerm(x, e, xi, T)[1];
			return ((1 - this->soil.Sr_w(x, e, xi, T)) * pow(p0/pC, lambda) + this->soil.Sr_n(x, e, xi, T));
		}

		double dSdP (double pC, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T) const
		{
			double lambda = this->soil.paramRelPerm(x, e, xi, T)[0];
			double p0 = this->soil.paramRelPerm(x, e, xi, T)[1];
			return (-(1 - this->soil.Sr_w(x, e, xi, T))*pow(p0/pC, lambda-1)*lambda*pow(1.0/pC, 2));
		}


		double krw (double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T) const
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


		double krn (double saturationN, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T) const
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

		std::vector<double> kr (const double saturationW, const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T) const
		{
			std::vector<double> kr(2);
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

		BrooksCoreyLaw(const Matrix2p<G,double>& s, bool lin = false)
			: RelPerm_pc<G>(s, false)
		{
		}

	private:
		static const double epsPC = 0.0001;
		static const double epsKr = 1e-15;
	};
}

#endif
