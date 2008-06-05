#ifndef DUNE_PROPERTIES_HH
#define DUNE_PROPERTIES_HH

/**
 * \ingroup material
 * \defgroup properties Fluid properties
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
  virtual double viscosity (double T=283.15, double p=1e5) const = 0;

  /** @brief density
   * @param T Temperature \f$ \left[ K \right] \f$
   * @param p Pressure \f$ \left[ Pa \right] \f$
   * @return density \f$ \left[ \frac{kg}{m^3} \right] \f$
   */
  virtual double density (double T=283.15, double p=1e5) const = 0;

  /** @brief residual saturation
   * @return residual saturation \f$ \left[ \, \right] \f$
   */
  virtual double Sr() const = 0;
  
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
  
  virtual ~Medium()
  { }
};

/** \ingroup properties
 * @brief Fluid properties of water
 */
class Water : public Medium
{
public:
	double viscosity (double T=283.15, double p=1e5) const
	{
		return 1e-3; //[kg/(ms)]
	}
	double density (double T=283.15, double p=1e5) const
	{
		return 1000.; // [kg/m^3]
	}
	double Sr() const
	{
		return 0.0;
	}
	double henry (double T=283.15) const
	{
		return (0.8942 + 1.47 * exp(-0.04394*T) )*1e-10; // [1/Pa]
	}
	double vaporPressure (double T=283.15) const
	{
		return 1228.; // [Pa]
	}
	double molarMass() const
	{
		return 0.018016; // [kg/mole]
	}
};

/** \ingroup properties
 * @brief Fluid properties of Air
 */
class Air : public Medium
{
public:
	double viscosity ( double T=283.15, double p=1e5) const
	{
		return 1.8e-5;//[kg/(ms)]
	}
	double density ( double T=283.15, double p=1e5) const
	{
		return 1.29; // [kg/m^3]
	}
	double Sr() const
	{
		return 0.1;
	}
	double molarMass() const
	{
		return 0.02896; // [kg/mole]
	}
};


class Brine : public Medium
{
public:
	double viscosity ( double T=283.15, double p=1e5) const
	{
	
		return 2.535e-4;//800e-3;//[kg/(ms)]
	}
	double density ( double T=283.15, double p=1e5) const
	{
		return 1045.0;//820.0; // [kg/m^3]
	}
	double Sr() const
	{
		return 0.0;
	}
	double molarMass() const
	{
		return 0; // [kg/mole]
	}
};

class Oil : public Medium
{
public:
	double viscosity ( double T=283.15, double p=1e5) const
	{
		return 1e-3;//800e-3;//[kg/(ms)]
	}
	double density ( double T=283.15, double p=1e5) const
	{
		return 1000.0;//820.0; // [kg/m^3]
	}
	double Sr() const
	{
		return 0.0;
	}
	double molarMass() const
	{
		return 0; // [kg/mole]
	}
};


/** \ingroup properties
 * @brief Uniform Fluid properties for testing and debugging.
 */
class Uniform : public Medium
{
public:
	double viscosity ( double T=283.15, double p=1e5) const
	{
		return 1.0;//[kg/(ms)]
	}
	double density ( double T=283.15, double p=1e5) const
	{
		return 1.0; // [kg/m^3]
	}
	double Sr() const
	{
		return 0.0;
	}
	double molarMass() const
	{
		return 1.0; // [kg/mole]
	}
};

/** \ingroup properties
 * @brief Fluid properties of DNAPL
 */
class DNAPL : public Medium
{
public:
	double viscosity ( double T=283.15, double p=1e5) const
	{
		return 5.7e-4;//[kg/(ms)]
	}
	double density ( double T=283.15, double p=1e5) const
	{
		return 1460.0; // [kg/m^3]
	}
	double Sr() const
	{
		return 0.0;
	}
	double molarMass() const
	{
		return 1.0; // [kg/mole]
	}
};

/** \ingroup properties
 * @brief Fluid properties of CO2 
 */
class CO2 : public Medium
{
public:
        double viscosity ( double T=432.0, double p=3.086e7) const // according to "webbook.nist.gov/chemistry/fluid/"
        {
                return 3.95e-5;//[kg/(ms)] // given in "CO2 leakage through an abandoned well.." A.Ebigbo,H.Class,R.Helmig 
        }
        double density ( double T=432.0, double p=3.086e7) const
        {
                return 479.0; // [kg/m^3] // given in "CO2 leakage through an abandoned well.." A.Ebigbo,H.Class,R.Helmig 
        }
        double Sr() const
        {
                return 0.0;
        }
        double molarMass() const
        {
                return 44.01; // [kg/mole]
        }
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

