#ifndef PHASEPROPERTIES_HH_
#define PHASEPROPERTIES_HH_

#include <dumux/material/constrel/constrelair.hh>
#include <dumux/material/constrel/constrelwater.hh>
#include <math.h>

/**
 * \ingroup properties
 * \author Jochen Fritz
 */

namespace Dune
{

/**\ingroup properties
 * @brief property class for gaseous phase of air/water mixture
 */
class G_air_water : public gas_gl
{
	virtual double density(double p, double T, double Xw) const // [kg / m^3]
	{
		double Rsm = R * (Xw / M_w + (1-Xw) / M_a); // medium specific gas constant
		return p / Rsm / T;
	}
	
	virtual double viscosity(double p, double T, double Xw) const // [kg / (m*s)]
	{
		double v_a = constRelAir.viscosity_air (T); // see constrelair.hh
		double v_w = constRelAir.visco_w_vap(T);    // see constrelair.hh
		FieldVector<double,2> X(Xw); X[1] = (1-Xw);
		X = X2x(X);
		X[0] *= sqrt(M_w); X[1] *= sqrt(M_a);
		return (v_w * X[0] + v_a * X[1]) / (X[0] + X[1]); // after Herning & Zipperer, 1936
	}
	
	virtual intEnergy(double p, double T, double Xw) const
	{
		return enthalpy(p,T,Xw) - p/density(p,T,Xw);
	}
	
	virtual enthalpy(double p, double T, double Xw) const
	{
		double H_a = 1.005 * (T - 273.15);
		double H_w = constRelAir.hsat(T);
		return Xw * H_w + (1-Xw) * H_a;
	}
	
	virtual double diffCoeff(double p, double T) const
	{
		// D ~ T^(3/2) / p see Atkins:Physical Chemistry p.778!
		// for H2O and O2: D(273.15 K, 1e5 Pa) = 2.25 e-5
		return 2.25e-5 * pow(T/273.15, 2/3) * 1e5 / p;
	}
	
	virtual double Xw_Max(double p, double T) const
	{
		double pwsat = constRelAir.pwsat(T);
		FieldVector<double,dim> x(min(pwsat / p, 1)); x[1] = 1-x[0];
		x = x2X(x);
		return x[0];
	}
private:
	ConstRelAir constRelAir;
	static const double M_w = 0.018016; //[kg / mole]
	static const double M_a = 0.02896; //[kg / mole]
	static const double R = 8.314472; // universal gas constant [J / (mole * K)]
};

} // namespace

#endif /*PHASEPROPERTIES_HH_*/
