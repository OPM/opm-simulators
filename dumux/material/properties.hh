#ifndef DUNE_PROPERTIES_HH
#define DUNE_PROPERTIES_HH
#include <dumux/material/constrel/constrelco2.hh>
#include <dumux/material/constrel/constrelwater.hh>
#include <dumux/material/constrel/constrelbrine.hh>
#include <dumux/material/constrel/constrelair.hh>
/**
 * \ingroup material
 * \defgroup properties Fluid and Soil properties
 */

/** \ingroup properties
 * @brief base class for fluid properties
 */
class Medium 
{
public:
  /** @brief kinematic viscosity
   * @param T Temperature \f$ \left[ K \right] \f$
   * @param p Pressure \f$ \left[ Pa \right] \f$
   * @return kinematic viscosity \f$ \left[ \frac{kg}{ms} \right] \f$
   */
  virtual double viscosity (double T=283.15, double p=1e5, double rho=0.) const = 0;

  /** @brief density
   * @param T Temperature \f$ \left[ K \right] \f$
   * @param p Pressure \f$ \left[ Pa \right] \f$
   * @return density \f$ \left[ \frac{kg}{m^3} \right] \f$
   */
  virtual double density (double T=283.15, double p=1e5, double X=1.) const = 0;

  /** @brief residual saturation
   * @return residual saturation \f$ \left[ \, \right] \f$
   */
  virtual double Sr() const = 0;

  /** @brief enthalpy
   * @param T Temperature \f$ \left[ K \right] \f$
   * @param p Pressure \f$ \left[ Pa \right] \f$
   * @return enthalpy \f$ \left[ \frac{J}{kg} \right] \f$
   */
  virtual double enthalpy (double T=283.15, double p=1e3, double X=1.) const
  {
	  return 0;
  }
  
  /**
   * @brief molar Mass of Component
   */
  virtual double molarMass() const = 0;
  
  /**
   * @brief vapor Pressure \f$ \left[ Pa \right] \f$
   * @param T Temperature \f$ \left[ K \right] \f$
   */
  double vaporPressure (double T) const
  {
	  return 0;
  }
  
  virtual double intEnergy(double T=283.15, double p=1e3) const
  {
	  return 0;
  }
 
  virtual ~Medium()
  { }
};

/** \ingroup properties
 * @brief base class for soil properties
 */
class Soil
{
public:
	double heatConductivity (double lDry = 0.32, double lSw = 2.7, double Sw = 1.)
	{
		double l_wurz;
		if(Sw<0.0) Sw = 0.0; /* ACHTUNG Regularisierung */
		if(Sw>1.) Sw = 1.; /* ACHTUNG Regularisierung */

		l_wurz = lDry + sqrt(Sw)*(lSw-lDry);

		if(isnan(l_wurz)) {
//			std::cout <<"isnan heatcondwurzel \n"<<std::endl;
			l_wurz = 0.0;
		}
		return(l_wurz);
	}
};

/** \ingroup properties
 * @brief Fluid properties of water
 */
class Water : public Medium
{
	ConstrelWater constRelWater;

public:
	Water(double Sr = 0.0, double constDensity = 0, 
			double constViscosity = 0, double constEnthalpy = 0)
	: Sr_(Sr), constDensity_(constDensity), constViscosity_(constViscosity), constEnthalpy_(constEnthalpy)
	{}
	
	double viscosity (double T=283.15, double p=1e5, double rho=0.) const
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
	double Sr() const
	{
		return Sr_;
	}
	double enthalpy (double T=283.15, double p=1e5, double X=1.) const
	{
		if (constEnthalpy_) 
			return constEnthalpy_;
		else {
			return constRelWater.enthalpy_water(T,p);
		}
	}
	double molarMass() const
	{
		return 0.018016; // [kg/mole]
	}
    double intEnergy( double T=283.15, double p=1e5) const
    {
    	double u;
    	double rho_mass = density(T,p);
    	double h = enthalpy(T,p);
     
    	u = h - (p / rho_mass);
    	return u;
    }
	double henry (double T=283.15) const
	{
		return (0.8942 + 1.47 * exp(-0.04394*T) )*1e-10; // [1/Pa]
	}
	double vaporPressure (double T=283.15) const
	{
		return 1228.; // [Pa]
	}

private:
   double Sr_;	
   double constDensity_;
   double constViscosity_;
   double constEnthalpy_;
};

/** \ingroup properties
 * @brief Fluid properties of Air
 */
class Air : public Medium
{
	ConstrelAir constRelAir;
	
public:
	Air(double Sr = 0.0, double constDensity = 0, 
			double constViscosity = 0, double constEnthalpy = 0)
	:Sr_(Sr), constDensity_(constDensity), constViscosity_(constViscosity), constEnthalpy_(constEnthalpy) 
	{}	
	double viscosity ( double T=283.15, double p=1e5, double rho=0.) const
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
			const double molarMassAir = molarMass();

			return constRelAir.rho_idGG_mass(T,p,molarMassAir); // [kg/m^3]
		}
	}
	double enthalpy (double T=283.15, double p=1e5, double Xwg=1.) const
	{
		if (constEnthalpy_) 
			return constEnthalpy_;
		else {
			return constRelAir.sp_enth2p2cni_g(T,p,Xwg);
		}
	}
	double Sr() const
	{
		return Sr_;
	}
	double molarMass() const
	{
		return 0.02896; // [kg/mole]
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
   double Sr_;	
   double constDensity_;
   double constViscosity_;
   double constEnthalpy_;
};


class Brine : public Medium
{
	ConstrelBrine constRelBrine;
	
public:
	Brine(double Sr = 0.0, double constDensity = 0, 
			double constViscosity = 0, double constEnthalpy = 0)
	: Sr_(Sr), constDensity_(constDensity), constViscosity_(constViscosity), constEnthalpy_(constEnthalpy) 
	{}
	
	double Salinity() const
	{
		return 0.1;
	}
	double viscosity ( double T=283.15, double p=1e5, double rho=0.) const
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
	double enthalpy (double T=283.15, double p=1e5, double X=1.) const
	{
		if (constEnthalpy_) 
			return constEnthalpy_;
		else {
			double S;
			S = Salinity();
			return constRelBrine.enthalpy_brine(T,p,S);
		}
	}
	double Sr() const
	{
		return Sr_;
	}
	double molarMass() const
	{
		return 0.018016; // [kg/mole]
	}
	double intEnergy(double T=283.15, double p=1e5) const
	{
		double intenergy;
		intenergy = enthalpy(T,p);
		return intenergy; 
	}

private:
   double Sr_;	
   double constDensity_;
   double constViscosity_;
   double constEnthalpy_;
};

class Oil : public Medium
{
public:
	Oil(double Sr = 0.0, double constDensity = 0, 
			double constViscosity = 0, double constEnthalpy = 0)
	: Sr_(Sr), constDensity_(constDensity), constViscosity_(constViscosity), constEnthalpy_(constEnthalpy) 
	{}
	
	double viscosity ( double T=283.15, double p=1e5, double rho=0.) const
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
	double enthalpy (double T=283.15, double p=1e5, double X=1.) const
	{
		if (constEnthalpy_) 
			return constEnthalpy_;
		else {
//			return constRelOil.enthalpy_brine(T,p,S);
		}
	}
	double Sr() const
	{
		return Sr_;
	}
	double molarMass() const
	{
		return 0; // [kg/mole]
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
   double Sr_;	
   double constDensity_;
   double constViscosity_;
   double constEnthalpy_;
};


/** \ingroup properties
 * @brief Uniform Fluid properties for testing and debugging.
 */
class Uniform : public Medium
{
public:
	Uniform(double Sr = 0.0, double constDensity = 0, 
			double constViscosity = 0, double constEnthalpy = 0)
	: Sr_(Sr), constDensity_(constDensity), constViscosity_(constViscosity), constEnthalpy_(constEnthalpy) 
	{}
	
	double viscosity ( double T=283.15, double p=1e5, double rho=0.) const
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
	double Sr() const
	{
		return Sr_;
	}
	double enthalpy (double T=283.15, double p=1e5, double X=1.) const
	{
		if (constEnthalpy_) 
			return constEnthalpy_;
		else {
			return 1.0;
		}
	}
	double molarMass() const
	{
		return 1.0; // [kg/mole]
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
   double Sr_;	
   double constDensity_;
   double constViscosity_;
   double constEnthalpy_;
};

/** \ingroup properties
 * @brief Fluid properties of DNAPL
 */
class DNAPL : public Medium
{
public:
	DNAPL(double Sr = 0.0, double constDensity = 0, 
			double constViscosity = 0, double constEnthalpy = 0)
	: Sr_(Sr), constDensity_(constDensity), constViscosity_(constViscosity), constEnthalpy_(constEnthalpy) 
	{}
	
	double viscosity ( double T=283.15, double p=1e5, double rho=0.) const
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
	double Sr() const
	{
		return Sr_;
	}
	double enthalpy (double T=283.15, double p=1e5, double X=1.) const
	{
		if (constEnthalpy_) 
			return constEnthalpy_;
		else {
//			return 1.0;
		}
	}
	double molarMass() const
	{
		return 1.0; // [kg/mole]
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
   double Sr_;	
   double constDensity_;
   double constViscosity_;
   double constEnthalpy_;
};

/** \ingroup properties
 * @brief Fluid properties of CO2 
 */
class CO2 : public Medium
{
 ConstrelCO2 constRelCO2;
	
 public:
		CO2(double Sr = 0.0, double constDensity = 0, 
				double constViscosity = 0, double constEnthalpy = 0)
		: Sr_(Sr), constDensity_(constDensity), constViscosity_(constViscosity), constEnthalpy_(constEnthalpy) 
		{}
		
		double density ( double T, double p, double X=1.) const
        {
			if (constDensity_) 
				return constDensity_;
			else 
				return constRelCO2.density(T,p);
		}
		double viscosity ( double T=432., double p=3.086e7, double rho=0.) const 
		{
			if (constViscosity_) 
				return constViscosity_;
			else 
				return constRelCO2.viscosity(T,p,rho);
		}
		double enthalpy ( double T=432., double p=3.086e7, double X=1.) const
		{
			if (constEnthalpy_) 
				return constEnthalpy_;
			else {
				return constRelCO2.enthalpy(T,p);
			}
		}
        double Sr() const
        {
            return Sr_;
        }
        double molarMass() const
        {
            return 0.04401; // [kg/mole]
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
    double Sr_;	
    double constDensity_;
    double constViscosity_;
    double constEnthalpy_;
};

/**\ingroup properties
 * @brief base class for the Henry coefficient of a binary gas/liquid system
 */
class Henry
{
public:	
	/** @brief Henry coefficient
	 * @param T Temperature \f$ \left[ K \right] \f$
	 * @return Henry coefficient \f$ \left[ 1/Pa \right] \f$
	 */
	virtual double operator() (double T) = 0; 
	
	virtual ~Henry() 
	{ }
};

/**\ingroup properties
 * @brief class for the Henry coefficient of a binary air/water system
 */
class HenryWaterAir : public Henry
{
public:
	virtual double operator() (double T)
	{
		return (0.8942 + 1.47 * exp(-0.04394*T) )*1e-10; // [1/Pa]
	}
};
#endif

