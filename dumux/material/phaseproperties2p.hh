// $Id$ 

#ifndef DUNE_PHASEPROPERTIES2P_HH
#define DUNE_PHASEPROPERTIES2P_HH

#include <dumux/material/property_baseclasses.hh>
#include <dumux/material/constrel/constrelco2.hh>
#include <dumux/material/constrel/constrelwater.hh>
#include <dumux/material/constrel/constrelbrine.hh>
#include <dumux/material/constrel/constrelair.hh>

namespace Dune
{
/** \ingroup properties
 * @brief Fluid properties of water
 */
class Water : public MediumNonIsothermal
{
	ConstrelWater constRelWater;

public:
	Water(double constDensity = 0,
			double constViscosity = 0, double constEnthalpy = 0)
	: constDensity_(constDensity), constViscosity_(constViscosity), constEnthalpy_(constEnthalpy)
	{}

	double viscosity (double T=283.15, double p=1e5, double X=1.) const
	{
		if (constViscosity_)
			return constViscosity_;
		else
			return constRelWater.viscosity_water(T,p); //[kg/(ms)]
	}

	double viscosity (double T=283.15, double p=1e5, double rho=0., double X=1.) const
	{
		if (constViscosity_)
			return constViscosity_;
		else
			return constRelWater.viscosity_water(T,p); //[kg/(ms)]
	}
	double density (double T=283.15, double p=1e5, double X=1.) const
	{
		if (constDensity_)
			return constDensity_;
		else
			return 1000.0; // assumed to be incompressible[kg/m^3]
	}

	double enthalpy (double T=283.15, double p=1e5) const
	{
		if (constEnthalpy_)
			return constEnthalpy_;
		else {
			return constRelWater.enthalpy_water(T,p);
		}
	}
    double intEnergy( double T=283.15, double p=1e5) const
    {
    	double u;
    	double rho_mass = density(T,p);
    	double h = enthalpy(T,p);

    	u = h - (p / rho_mass);
    	return u;
    }

private:
   double constDensity_;
   double constViscosity_;
   double constEnthalpy_;
};

/** \ingroup properties
 * @brief Fluid properties of Air
 */
class Air : public MediumNonIsothermal
{
	ConstrelAir constRelAir;

public:
	Air(double constDensity = 0,
			double constViscosity = 0, double constEnthalpy = 0)
	:constDensity_(constDensity), constViscosity_(constViscosity), constEnthalpy_(constEnthalpy)
	{}

	double viscosity ( double T=283.15, double p=1e5, double X=1.) const
	{
		if (constViscosity_)
			return constViscosity_;
		else
			return constRelAir.viscosity_air(T); //[kg/(ms)]
	}

	double viscosity ( double T=283.15, double p=1e5, double rho=0., double X=1.) const
	{
		if (constViscosity_)
			return constViscosity_;
		else
			return constRelAir.viscosity_air(T); //[kg/(ms)]
	}
	double density ( double T=283.15, double p=1e5, double X=1.) const
	{
		if (constDensity_)
			return constDensity_;
		else {
			const double molarMassAir = 0.02896; // [kg/mole]

			return constRelAir.rho_idGG_mass(T,p,molarMassAir); // [kg/m^3]
		}
	}
	double enthalpy (double T=283.15, double p=1e5) const
	{
		if (constEnthalpy_)
			return constEnthalpy_;
		else {
			return constRelAir.sp_enth2p2cni_g(T,p,1);
		}
	}

    double intEnergy( double T=283.15, double p=1e5) const
    {
    	double u;
    	double rho_mass = density(T,p);
    	double h = enthalpy(T,p);

    	u = h - (p / rho_mass);
    	return u;
    }

private:
   double constDensity_;
   double constViscosity_;
   double constEnthalpy_;
};


class Brine : public MediumNonIsothermal
{
	ConstrelBrine constRelBrine;

public:
	Brine(double constDensity = 0,
			double constViscosity = 0, double constEnthalpy = 0)
	:constDensity_(constDensity), constViscosity_(constViscosity), constEnthalpy_(constEnthalpy)
	{}

	double viscosity ( double T=283.15, double p=1e5, double X=1.) const
	{
		if (constViscosity_)
			return constViscosity_;
		else {
			double S;
			S = Salinity();
			return constRelBrine.viscosity_brine(T,S);
		}

//  		 return 2.535e-4; // [kg/(ms)]
	}

	double Salinity() const
	{
		return 0.1;
	}
	double viscosity ( double T=283.15, double p=1e5, double rho=0., double X=1.) const
	{
		if (constViscosity_)
			return constViscosity_;
		else {
			double S;
			S = Salinity();
			return constRelBrine.viscosity_brine(T,S);
		}

//  		 return 2.535e-4; // [kg/(ms)]
	}
	double density ( double T=283.15, double p=1e5, double X=1.) const
	{
		if (constDensity_)
			return constDensity_;
		else {
			double S, x_CO2_w;
			x_CO2_w = 0.0;
			S = Salinity();
			return constRelBrine.mass_density_brine_CO2(T,p,S,x_CO2_w);
		}
	}
	double enthalpy (double T=283.15, double p=1e5) const
	{
		if (constEnthalpy_)
			return constEnthalpy_;
		else {
			double S;
			S = Salinity();
			return constRelBrine.enthalpy_brine(T,p,S);
		}
	}

	double intEnergy(double T=283.15, double p=1e5) const
	{
		double intenergy;
		intenergy = enthalpy(T,p);
		return intenergy;
	}

private:

   double constDensity_;
   double constViscosity_;
   double constEnthalpy_;
};

class Oil : public MediumNonIsothermal
{
public:
	Oil(double constDensity = 0,
			double constViscosity = 0, double constEnthalpy = 0)
	: constDensity_(constDensity), constViscosity_(constViscosity), constEnthalpy_(constEnthalpy)
	{}

	double viscosity ( double T=283.15, double p=1e5, double X=1.) const
	{
		if (constViscosity_)
			return constViscosity_;
		else
			return 1e-3;//800e-3;//[kg/(ms)]
	}

	double viscosity ( double T=283.15, double p=1e5, double rho=0., double X=1.) const
	{
		if (constViscosity_)
			return constViscosity_;
		else
			return 1e-3;//800e-3;//[kg/(ms)]
	}
	double density ( double T=283.15, double p=1e5, double X=1.) const
	{
		if (constDensity_)
			return constDensity_;
		else
			return 1000.0;//820.0; // [kg/m^3]
	}
	double enthalpy (double T=283.15, double p=1e5) const
	{
		if (constEnthalpy_)
			return constEnthalpy_;
		else {
//			return constRelOil.enthalpy_brine(T,p,S);
		}
	}

    double intEnergy( double T=283.15, double p=1e5) const
    {
    	double u;
    	double rho_mass = density(T,p);
    	double h = enthalpy(T,p);

    	u = h - (p / rho_mass);
    	return u;
    }

private:
   double constDensity_;
   double constViscosity_;
   double constEnthalpy_;
};


/** \ingroup properties
 * @brief Uniform Fluid properties for testing and debugging.
 */
class Uniform : public MediumNonIsothermal
{
public:
	Uniform(double constDensity = 0,
			double constViscosity = 0, double constEnthalpy = 0)
	: constDensity_(constDensity), constViscosity_(constViscosity), constEnthalpy_(constEnthalpy)
	{}

	double viscosity ( double T=283.15, double p=1e5, double X=1.) const
	{
		if (constViscosity_)
			return constViscosity_;
		else
			return 1.0;//[kg/(ms)]
	}

	double viscosity ( double T=283.15, double p=1e5, double rho=0., double X=1.) const
	{
		if (constViscosity_)
			return constViscosity_;
		else
			return 1.0;//[kg/(ms)]
	}
	double density ( double T=283.15, double p=1e5, double X=1.) const
	{
		if (constDensity_)
			return constDensity_;
		else
			return 1.0; // [kg/m^3]
	}

	double enthalpy (double T=283.15, double p=1e5) const
	{
		if (constEnthalpy_)
			return constEnthalpy_;
		else {
			return 1.0;
		}
	}
    double intEnergy( double T=283.15, double p=1e5) const
    {
    	double u;
    	double rho_mass = density(T,p);
    	double h = enthalpy(T,p);

    	u = h - (p / rho_mass);
    	return u;
    }

private:
   double constDensity_;
   double constViscosity_;
   double constEnthalpy_;
};

/** \ingroup properties
 * @brief Fluid properties of DNAPL
 */
class DNAPL : public MediumNonIsothermal
{
public:
	DNAPL(double constDensity = 0,
			double constViscosity = 0, double constEnthalpy = 0)
	: constDensity_(constDensity), constViscosity_(constViscosity), constEnthalpy_(constEnthalpy)
	{}

	double viscosity ( double T=283.15, double p=1e5, double X=1.) const
	{
		if (constViscosity_)
			return constViscosity_;
		else
			return 5.7e-4;//[kg/(ms)]
	}

	double viscosity ( double T=283.15, double p=1e5, double rho=0., double X=1.) const
	{
		if (constViscosity_)
			return constViscosity_;
		else
			return 5.7e-4;//[kg/(ms)]
	}
	double density ( double T=283.15, double p=1e5, double X=1.) const
	{
		if (constDensity_)
			return constDensity_;
		else
			return 1460.0; // [kg/m^3]
	}
	double enthalpy (double T=283.15, double p=1e5) const
	{
		if (constEnthalpy_)
			return constEnthalpy_;
		else {
//			return 1.0;
		}
	}
    double intEnergy( double T=283.15, double p=1e5) const
    {
    	double u;
    	double rho_mass = density(T,p);
    	double h = enthalpy(T,p);

    	u = h - (p / rho_mass);
    	return u;
    }

private:
   double constDensity_;
   double constViscosity_;
   double constEnthalpy_;
};

/** \ingroup properties
 * @brief Fluid properties of CO2
 */
class CO2 : public MediumNonIsothermal
{
 ConstrelCO2 constRelCO2;

 public:
		CO2(double constDensity = 0,
				double constViscosity = 0, double constEnthalpy = 0)
		: constDensity_(constDensity), constViscosity_(constViscosity), constEnthalpy_(constEnthalpy)
		{}

		double viscosity ( double T=432., double p=3.086e7, double X=1.) const
		{
			if (constViscosity_)
				return constViscosity_;
			else
				return 0;
		}

		double density ( double T, double p, double X=1.) const
        {
			if (constDensity_)
				return constDensity_;
			else
				return constRelCO2.density(T,p);
		}
		double viscosity ( double T=432., double p=3.086e7, double rho=0., double X=1.) const
		{
			if (constViscosity_)
				return constViscosity_;
			else
				return constRelCO2.viscosity(T,p,rho);
		}
		double enthalpy ( double T=432., double p=3.086e7) const
		{
			if (constEnthalpy_)
				return constEnthalpy_;
			else {
				return constRelCO2.enthalpy(T,p);
			}
		}

        double intEnergy( double T=432, double p=3.086e7) const
        {
        	double u;
        	double rho_mass = density(T,p);
        	double h = enthalpy(T,p);

        	u = h - (p / rho_mass);
        	return u;
        }

 private:
    double constDensity_;
    double constViscosity_;
    double constEnthalpy_;
};
}
#endif

