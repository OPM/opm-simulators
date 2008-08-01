#ifndef PHASES_HH_
#define PHASES_HH_

#include <dune/common/fvector.hh>

/**
 * \ingroup material
 * \defgroup phase properties
 * \author Jochen Fritz
 */

namespace Dune
{
/** @brief Base class for liquid phase of a binary gas liquid mixture
 * Base class for a liquid phase made up of two components, where one
 * is a solved gas.
 * The components are denoted by 
 * w for the main (liquid) component and
 * a for the dissolved (gaseous) component.
 */
class 2liquid_gl
{
	//! Phase density in [kg / m^3]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 *  @param Xa mass fraction of dissolved component
	 */
	virtual double density(double p, double T, double Xa) = 0; // [kg / m^3]

	//! Dynamic viscosity in [kg / (m*s)]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 *  @param Xa mass fraction of dissolved component
	 */
	virtual double viscosity(double p, double T, double Xa) = 0; // [kg / (m*s)]
	
	//! Dynamic viscosity in [kg / (m*s)]
	/** by specification of the phase density, the internal density calculation can be avoided
	 *  @param p pressure [Pa]
	 *  @param T temperature [K]
	 *  @param rho Phase density [kg / m^3]
	 *  @param Xa mass fraction of dissolved component
	 */
	virtual double viscosity(double p, double T, double rho, double Xa) = 0; // [kg / (m*s)]
	
	//! Internal energy in [N * m]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 *  @param Xa mass fraction of dissolved component
	 */
	virtual double intEnergy(double p, double T, double Xa) = 0; // [N * m]
	
	//! Diffusion coefficient for component a [m^2 / s]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 */
	virtual double diffCoeff(double p, double T) = 0; // [m^2 / s]
	
	//! Solubility of component a [kg / kg]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 */
	virtual double Xa_Max(double p, double T) = 0; // [kg / kg]
	
	//! vapour pressure of liquid component in [Pa]
	/** @param T temperature in [K]
	*/
	virtual double p_vap(double T) = 0; // [Pa]
	
	//! Henry coefficient of component a [1 / Pa]
	/** @param T temperature [K]
	 */
	virtual double henry(double T) = 0; // [1 / Pa]
	
	//! Boiling point temperature of liquid components in [K]
	/** @param pressure in [Pa] 
	 */
	virtual double T_vap(double p) = 0; // [K]
	
	//! Molar mass of liquid component
	double molarMass_w()
	{
		return M_w;
	}
	
	//! Molar mass of gaseous solute
	double molarMass_a()
	{
		return M_a;
	}
	
	//! Molar mass of both components where the first entry is the liquid component
	FieldVector<double,2> molarMass()
	{
		FieldVector<double,2> M;
		M[0] = M_w;
		M[1] = M_a;
		return M;
	}
	
	//! Conversion from mole fractions to mass fractions
	/** @param x mole fractions of the components. The first entry represents the liquid component
	 */
	FieldVector<double,2> x2X(FieldVector<double,2> x)
	{
		if (x[0]+x[1] != 1.) DUNE_THROw(MathError, "mole fractions do not sum up to unity!");
		FieldVector<double,2> X;
		X[0] = 1 / (1 + x[1] * M_a / (x[0] * M_w));
		X[1] = 1 - X[0];
		return X;
	}
	
	//! Conversion from mass fractions to mole fractions
	/** @param X mass fractions of the components. The first entry represents the liquid component
	 */
	FieldVector<double,2> X2x(FieldVector<double,2> X)
	{
		if (X[0]+X[1] != 1.) DUNE_THROw(MathError, "mole fractions do not sum up to unity!");
		FieldVector<double,2> X;
		x[0] = 1 / (1 + X[1] / M_a / (X[0] / M_w));
		x[1] = 1 - x[0];
		return x;
	}
	
private:
	static const double M_w;
	static const double M_a;
};


/** @brief Base class for the liquid phase of a gas liquid mixture
 * Base class for a liquid phase made up of two sorts of components:
 * The liquid components, denoted by "w" and the dissolved
 * gaseous components, denoted by "a"
 * Template parameters:
 * n_w : number of liquid components
 * n_a : number of gaseous components
 */
template <int n_w, int n_a>
class nliquid_gl
{
	//! Phase density in [kg / m^3]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 *  @param Xa mass fraction of dissolved component
	 */
	virtual double density(double p, double T, FieldVector<double,(n_w + n_a)> X) = 0;
	
	//! Dynamic viscosity in [kg / (m*s)]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 *  @param Xa mass fraction of dissolved component
	 */
	virtual double viscosity(double p, double T, FieldVector<double,(n_w + n_a)>) = 0;
	
	//! Dynamic viscosity in [kg / (m*s)]
	/** by specification of the phase density, the internal density calculation can be avoided
	 *  @param p pressure [Pa]
	 *  @param T temperature [K]
	 *  @param rho Phase density [kg / m^3]
	 *  @param Xa mass fraction of dissolved component
	 */	
	virtual double viscosity(double p, double T, FieldVector<double,(n_w + n_a)>, double rho) = 0;
	
	//! Internal energy in [N * m]
	/** @param p pressure [Pa]
	 *  @param T temperature [K]
	 *  @param Xa mass fraction of dissolved component
	 */
	virtual double intEnergy(double p, double T, FieldVector<double,(n_w + n_a)>) = 0;

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
		return M_
	}
	
	//! Conversion from mole fractions to mass fractions
	/** @param x mole fractions of components (have to be in same order as molar masses in private vector M_ !)
	 */
	FieldVector<double,(n_w+n_a)> x2X(FieldVector<double,(n_w+n_a) x)
	{
		FieldVector<double,(n_w+n_a)> X(x);
		int size = n_w + n_a;
		double divisor;
		for (int i = 0; i < size; i++)
		{
			X[i] *= M[i];
			divisor += X[i];
		}
		X[i] /= divisor;
		return X;
	}
	
	//! Conversion from mass fractions to mole fractions
	/** @param X mass fractions of components (have to be in same order as molar masses in private vector M_ !)
	 */
	FieldVector<double,(n_w+n_a)> X2x(FieldVector<double,(n_w+n_a) X)
	{
		FieldVector<double,(n_w+n_a)> x(X);
		int size = n_w + n_a;
		double divisor;
		for (int i = 0; i < size; i++)
		{
			x[i] /= M[i];
			divisor += x[i];
		}
		x[i] /= divisor;
		return x;
	}
	
private:
	FieldVector<double,(n_w + n_a)> M_
};
}
#endif /*PHASES_HH_*/
