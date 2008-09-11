// $Id$

#ifndef PROPERTY_BASECLASSES
#define PROPERTY_BASECLASSES

#include <dune/common/fvector.hh>
#include <vector>

/**
 * \ingroup properties
 * \author Jochen Fritz
 */

namespace Dune
{
/**\ingroup properties
 * @brief Base class for matrix properties in case of two fluid phases
 * phases are denoted by
 * w for wetting phase
 * n for non-wetting phase
 */
template<class G, class RT>
class Matrix2p
{
public:
typedef	typename G::ctype DT;
	enum
	{	n=G::dimension, m=2};
	typedef typename G::Traits::template Codim<0>::Entity Entity;

	//! provides flags for the relative permeability and capillary pressure model
	enum modelFlag
	{
		linear = 0,
		brooks_corey = 1,
		van_genuchten = 2,
		auxiliary1 = 3,
		auxiliary2 = 4,
		auxiliary3 = 5,
	};

	/** @brief Permeability tensor
	 * @param x position in global coordinates
	 * @param e codim 0 entity for which the value is sought
	 * @param xi position in local coordinates in e
	 */
	virtual FieldMatrix<DT,n,n> K (const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) = 0;

	/**@brief matrix porosity
	 * @param x position in global coordinates
	 * @param e codim 0 entity for which the value is sought
	 * @param xi position in local coordinates in e
	 */
	virtual double porosity(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const = 0;

	/**@brief Wetting phase residual saturation
	 * @param x position in global coordinates
	 * @param e codim 0 entity for which the value is sought
	 * @param xi position in local coordinates in e
	 * @param T temperature in [K]
	 */
	virtual double Sr_w(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T = 283.15) const = 0;

	/**@brief Nonwetting phase residual saturation
	 * @param x position in global coordinates
	 * @param e codim 0 entity for which the value is sought
	 * @param xi position in local coordinates in e
	 * @param T temperature in [K]
	 */
	virtual double Sr_n(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T = 283.15) const = 0;

	/**@brief Matrix heat capacity in [kJ / (m^3 K)]
	 * ATTENTION: define heat capacity per cubic meter! Be sure, that it corresponds to porosity!
	 * Best thing will be to define heatCap = (specific heatCapacity of material) * density * porosity
	 * @param x position in global coordinates
	 * @param e codim 0 entity for which the value is sought
	 * @param xi position in local coordinates in e
	 */
	/* ATTENTION: define heat capacity per cubic meter! Be sure, that it corresponds to porosity!
	 * Best thing will be to define heatCap = (specific heatCapacity of material) * density * porosity*/
	virtual double heatCap(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
	{
		DUNE_THROW(NotImplemented, "heat capacity function not implemented!");
	}

	/**@brief Heat conductivity of matrix AND fluid phases [ W / (m * K)]
	 * @param x position in global coordinates
	 * @param e codim 0 entity for which the value is sought
	 * @param xi position in local coordinates in e
	 * @param sat wetting Phase saturation
	 */
	virtual double heatCond(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double sat) const
	{
		DUNE_THROW(NotImplemented, "heat conductivity function not implemented!");
	}

	/**@brief parameters for relative permeabilty models
	 * @param x position in global coordinates
	 * @param e codim 0 entity for which the value is sought
	 * @param xi position in local coordinates in e
	 * @param T Temperature
	 */
	virtual std::vector<double> paramRelPerm(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T = 283.15) const = 0;

	/**@brief Flag for determining the relative permeability model
	 * @param x position in global coordinates
	 * @param e codim 0 entity for which the value is sought
	 * @param xi position in local coordinates in e
	 */
	virtual modelFlag relPermFlag(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
	{
		return linear;
	}

	virtual ~Matrix2p()
	{}
};

/**ingroup properties
 * @brief baseclass for fluid phases.
 */
class Medium
{
public:
	/** @brief dynamic viscosity in [kg / (m*s)]
	 * @param T Temperature \f$ \left[ K \right] \f$
	 * @param p Pressure \f$ \left[ Pa \right] \f$
	 *  @param Xa mass fraction of dissolved component
	 * @return kinematic viscosity \f$ \left[ \frac{kg}{ms} \right] \f$
	 */
	virtual double viscosity (double T = 283.15, double p = 1e5, double X = 0.) const = 0;

	//! dynamic viscosity in [kg / (m*s)]
	/** by specification of the phase density, the internal density calculation can be avoided
	 *  @param p pressure [Pa]
	 *  @param T temperature [K]
	 *  @param rho Phase density [kg / m^3]
	 *  @param X mass fraction of dissolved component
	 */
	virtual double viscosityCO2(double T, double p, double rho, double X = 0.) const // [kg / (m*s)]
	{
		return 0;
	}

	/** @brief density
	 * @param T Temperature \f$ \left[ K \right] \f$
	 * @param p Pressure \f$ \left[ Pa \right] \f$
	 * @param X mass fraction of dissolved component
	 * @return density \f$ \left[ \frac{kg}{m^3} \right] \f$
	 */
	virtual double density (double T = 283.15, double p = 1e5, double X = 0.) const = 0;


	virtual ~Medium()
	{}
};

/** \ingroup properties
 * @brief base class for fluid properties including properties needed for the calculation of  non-isothermal processes
 */
class MediumNonIsothermal : public Medium
{
public:
	/** @brief enthalpy
	 * @param T Temperature \f$ \left[ K \right] \f$
	 * @param p Pressure \f$ \left[ Pa \right] \f$
	 * @return enthalpy \f$ \left[ \frac{J}{kg} \right] \f$
	 */
	virtual double enthalpy (double T=283.15, double p=1e5, double X = 1.) const = 0;

	virtual double intEnergy(double T=283.15, double p=1e5, double X = 1.) const =0;

	virtual ~MediumNonIsothermal()
	{}
};

/**\ingroup properties
 * @brief Base class for liquid phase of a binary gas liquid mixture
 * Base class for a liquid phase made up of two components, where one
 * is a solved gas.
 * The components are denoted by
 * w for the main (liquid) component and
 * a for the dissolved (gaseous) component.
 */
class Liquid_GL : public MediumNonIsothermal
{
public:

	/** @brief dynamic viscosity in [kg / (m*s)]
	 * @param T Temperature \f$ \left[ K \right] \f$
	 * @param p Pressure \f$ \left[ Pa \right] \f$
	 *  @param Xa mass fraction of dissolved component
	 * @return kinematic viscosity \f$ \left[ \frac{kg}{ms} \right] \f$
	 */
	virtual double viscosity (double T, double p, double X = 0.) const = 0;

	/** @brief density
	 * @param T Temperature \f$ \left[ K \right] \f$
	 * @param p Pressure \f$ \left[ Pa \right] \f$
	 * @param X mass fraction of dissolved component
	 * @return density \f$ \left[ \frac{kg}{m^3} \right] \f$
	 */
	virtual double density (double T, double p, double X = 0.) const = 0;

	/** @brief enthalpy
	 * @param T Temperature \f$ \left[ K \right] \f$
	 * @param p Pressure \f$ \left[ Pa \right] \f$
	 * @return enthalpy \f$ \left[ \frac{J}{kg} \right] \f$
	 */
	virtual double enthalpy (double T, double p, double Xa = 0.) const = 0;

	//! Specific internal energy in [N * m / kg]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 *  @param Xa mass fraction of dissolved component
	 */
	virtual double intEnergy(double T, double p, double Xa = 0.) const = 0; // [N * m / kg]

	//! Diffusion coefficient for component a [m^2 / s]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 */
	virtual double diffCoeff(double T, double p) const = 0; // [m^2 / s]

	//! Solubility of component a [kg / kg]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 */
	virtual double Xa_Max(double T, double p) const = 0; // [kg / kg]

	//! vapour pressure of liquid component in [Pa]
	/** @param T temperature in [K]
	 */
	virtual double p_vap(double T) const = 0; // [Pa]

	//! Henry coefficient of component a [1 / Pa]
	/** @param T temperature [K]
	 */
	virtual double henry(double T) const = 0; // [1 / Pa]

	//! Boiling point temperature of liquid components in [K]
	/** @param pressure in [Pa]
	 */
	virtual double T_vap(double p) const = 0; // [K]

	//! Molar mass of liquid component
	virtual double molarMass_w()
	{
		return M_w;
	}

	//! Molar mass of gaseous solute
	virtual double molarMass_a()
	{
		return M_a;
	}

	//! Molar mass of both components where the first entry is the liquid component
	virtual FieldVector<double,2> molarMass()
	{
		FieldVector<double,2> M;
		M[0] = M_w;
		M[1] = M_a;
		return M;
	}

	//! Conversion from mole fractions to mass fractions
	/** @param x mole fractions of the components. The first entry represents the liquid component
	 */
	virtual FieldVector<double,2> x2X(FieldVector<double,2> x) const
	{
		if (x[0]+x[1] != 1.) DUNE_THROW(Dune::MathError, "mole fractions do not sum up to unity!");
		FieldVector<double,2> X;
		X[0] = 1 / (1 + x[1] * M_a / (x[0] * M_w));
		X[1] = 1 - X[0];
		return X;
	}

	//! Conversion from mass fractions to mole fractions
	/** @param X mass fractions of the components. The first entry represents the liquid component
	 */
	virtual FieldVector<double,2> X2x(FieldVector<double,2> X) const
	{
		if (X[0]+X[1] != 1.) DUNE_THROW(Dune::MathError, "mole fractions do not sum up to unity!");
		FieldVector<double,2> x;
		x[0] = 1 / (1 + X[1] / M_a / (X[0] / M_w));
		x[1] = 1 - x[0];
		return x;
	}

	Liquid_GL() : M_w(1.0), M_a(1.0)
	{
	}

protected:
	double M_w;
	double M_a;
};

/**\ingroup properties
 * @brief Base class for the liquid phase of a gas liquid mixture
 * Base class for a liquid phase made up of two sorts of components:
 * The liquid components, denoted by "w" and the dissolved
 * gaseous components, denoted by "a"
 * Template parameters:
 * n_w : number of liquid components
 * n_a : number of gaseous components
 */
template <int n_w, int n_a>
class NLiquid_GL
{
	//! Phase density in [kg / m^3]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 *  @param X phase composition: mass fractions of components
	 */
	virtual double density(double T, double p, FieldVector<double,(n_w + n_a)> X) = 0;

	//! Dynamic viscosity in [kg / (m*s)]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 *  @param X phase composition: mass fractions of components
	 */
	virtual double viscosity(double T, double p, FieldVector<double,(n_w + n_a)> X) = 0;

	/** @brief enthalpy
	 * @param T Temperature \f$ \left[ K \right] \f$
	 * @param p Pressure \f$ \left[ Pa \right] \f$
	 * @return enthalpy \f$ \left[ \frac{J}{kg} \right] \f$
	 */
	virtual double enthalpy (double T, double p, double Xa) const = 0;

	//! Specific internal energy in [N * m / kg]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 *  @param X phase composition: mass fractions of components
	 */
	virtual double intEnergy(double T, double p, FieldVector<double,(n_w + n_a)> X) = 0;

	//! Diffusion coefficient for component a [m^2 / s]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 */
	virtual FieldVector<double,(n_w + n_a)> diffCoeff(double p, double T) = 0;

	//! Solubility of component a [kg / kg]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 */
	virtual FieldVector<double,n_a> Xa_Max(double p, double T, FieldVector<double,n_w> Xw) = 0;

	//! vapour pressure of liquid components in [Pa]
	/** @param T temperature in [K]
	 */
	virtual FieldVector<double, n_w> p_vap(double T) = 0; // [Pa]

	//! Boiling point temperature of liquid components in [K]
	/** @param pressure in [Pa]
	 */
	virtual FieldVector<double, n_w> T_vap(double p) = 0; // [K]

	//! Henry coefficient of component a [1 / Pa]
	/** @param T temperature [K]
	 */
	virtual FieldVector<double,n_a> henry(double T, FieldVector<double,n_w> Xw) = 0;

	//! Molar masses of components
	FieldVector<double,(n_w+n_a)> molarMass()
	{
		return M_;
	}

	//! Conversion from mole fractions to mass fractions
	/** @param x mole fractions of components (have to be in same order as molar masses in private vector M_ !)
	 */
	FieldVector<double,(n_w+n_a)> x2X(FieldVector<double,(n_w+n_a)> x)
	{
		FieldVector<double,(n_w+n_a)> X(x);
		int size = n_w + n_a;
		double divisor;
		for (int i = 0; i < size; i++)
		{
			X[i] *= M_[i];
			divisor += X[i];
		}
		X /= divisor;
		return X;
	}

	//! Conversion from mass fractions to mole fractions
	/** @param X mass fractions of components (have to be in same order as molar masses in private vector M_ !)
	 */
	FieldVector<double,(n_w+n_a)> X2x(FieldVector<double,(n_w+n_a)> X)
	{
		FieldVector<double,(n_w+n_a)> x(X);
		int size = n_w + n_a;
		double divisor;
		for (int i = 0; i < size; i++)
		{
			x[i] /= M_[i];
			divisor += x[i];
		}
		x /= divisor;
		return x;
	}

private:
	Dune::FieldVector<double,(n_w + n_a)> M_;
};

/**\ingroup properties
 * @brief Base class for the gas phase of a binary gas liquid mixture
 * Base class for a gaseous phase made up of two components, where one
 * is a vaporized liquid.
 * The components are denoted by
 * w for the dissolved (liquid) component and
 * a for the main (gaseous) component.
 */
class Gas_GL : public MediumNonIsothermal
{
public:
	/** @brief dynamic viscosity in [kg / (m*s)]
	 * @param T Temperature \f$ \left[ K \right] \f$
	 * @param p Pressure \f$ \left[ Pa \right] \f$
	 *  @param Xa mass fraction of dissolved component
	 * @return kinematic viscosity \f$ \left[ \frac{kg}{ms} \right] \f$
	 */
	virtual double viscosity (double T, double p, double X = 0.) const = 0;

	/** @brief density
	 * @param T Temperature \f$ \left[ K \right] \f$
	 * @param p Pressure \f$ \left[ Pa \right] \f$
	 * @param X mass fraction of dissolved component
	 * @return density \f$ \left[ \frac{kg}{m^3} \right] \f$
	 */
	virtual double density (double T, double p, double X = 0.) const = 0;

	//! Specific internal energy in [N * m / kg]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 *  @param Xw mass fraction of dissolved component
	 */
	virtual double intEnergy(double T, double p, double Xw = 0.) const = 0; // [N * m / kg]

	//! Specific enthalpy in [N * m / kg]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 *  @param Xw mass fraction of dissolved component
	 */
	virtual double enthalpy(double T, double p, double Xw = 0.) const = 0; // [N * m / kg]

	//! Diffusion coefficient for component a [m^2 / s]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 */
	virtual double diffCoeff(double T, double p) const = 0; // [m^2 / s]

	//! Solubility of component w [kg / kg]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 */
	virtual double Xw_Max(double T, double p) const = 0; // [kg / kg]

	//! Molar mass of liquid solute
	double molarMass_w() const
	{
		return M_w;
	}

	//! Molar mass of gaseous component
	double molarMass_a() const
	{
		return M_a;
	}

	//! Molar mass of both components where the first entry is the liquid component
	FieldVector<double,2> molarMass() const
	{
		FieldVector<double,2> M;
		M[0] = M_w;
		M[1] = M_a;
		return M;
	}

	//! Conversion from mole fractions to mass fractions
	/** @param x mole fractions of the components. The first entry represents the liquid component
	 */
	FieldVector<double,2> x2X(FieldVector<double,2> x) const
	{
		if (x[0]+x[1] != 1.) DUNE_THROW(MathError, "mole fractions do not sum up to unity!");
		FieldVector<double,2> X;
		X[0] = 1 / (1 + x[1] * M_a / (x[0] * M_w));
		X[1] = 1 - X[0];
		return X;
	}

	//! Conversion from mass fractions to mole fractions
	/** @param X mass fractions of the components. The first entry represents the liquid component
	 */
	FieldVector<double,2> X2x(FieldVector<double,2> X) const
	{
		if (X[0]+X[1] != 1.) DUNE_THROW(MathError, "mass fractions do not sum up to unity!");
		FieldVector<double,2> x;
		x[0] = 1 / (1 + X[1] / M_a / (X[0] / M_w));
		x[1] = 1 - x[0];
		return x;
	}

private:
	static const double M_w = 1;
	static const double M_a = 1;
};

/**\ingroup properties
 * @brief Base class for the gaseous phase of a gas liquid mixture
 * Base class for a gaseous phase made up of two sorts of components:
 * The gaseous components, denoted by "a" and the
 * dissolved liquid components, denoted by "w"
 * Template parameters:
 * n_w : number of liquid components
 * n_a : number of gaseous components
 */
template <int n_w, int n_a>
class NGas_GL
{
	//! Phase density in [kg / m^3]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 *  @param X phase composition: mass fractions of components
	 */
	virtual double density(double T, double p, FieldVector<double,(n_w + n_a)> X) = 0;

	//! Dynamic viscosity in [kg / (m*s)]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 *  @param X phase composition: mass fractions of components
	 */
	virtual double viscosity(double T, double p, FieldVector<double,(n_w + n_a)> X) = 0;

	//! Specific internal energy in [N * m / kg]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 *  @param X phase composition: mass fractions of components
	 */
	virtual double intEnergy(double T, double p, FieldVector<double,(n_w + n_a)> X) = 0;

	//! Diffusion coefficient for component a [m^2 / s]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 */
	virtual FieldVector<double,(n_w + n_a)> diffCoeff(double p, double T) = 0;

	//! Solubility of component a [kg / kg]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 */
	virtual FieldVector<double,n_w> Xw_Max(double p, double T, FieldVector<double,n_w> Xw) = 0;

	//! Molar masses of components
	FieldVector<double,(n_w+n_a)> molarMass()
	{
		return M_;
	}

	//! Conversion from mole fractions to mass fractions
	/** @param x mole fractions of components (have to be in same order as molar masses in private vector M_ !)
	 */
	FieldVector<double,(n_w+n_a)> x2X(FieldVector<double,(n_w+n_a)> x)
	{
		FieldVector<double,(n_w+n_a)> X(x);
		int size = n_w + n_a;
		double divisor;
		for (int i = 0; i < size; i++)
		{
			X[i] *= M_[i];
			divisor += X[i];
		}
		X /= divisor;
		return X;
	}

	//! Conversion from mass fractions to mole fractions
	/** @param X mass fractions of components (have to be in same order as molar masses in private vector M_ !)
	 */
	FieldVector<double,(n_w+n_a)> X2x(FieldVector<double,(n_w+n_a)> X)
	{
		FieldVector<double,(n_w+n_a)> x(X);
		int size = n_w + n_a;
		double divisor;
		for (int i = 0; i < size; i++)
		{
			x[i] /= M_[i];
			divisor += x[i];
		}
		x /= divisor;
		return x;
	}

private:
	FieldVector<double,(n_w + n_a)> M_;
};
} // end namespace
#endif /*PROPERTY_BASECLASSES*/
