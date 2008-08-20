#ifndef MATRIXPROPERTIES
#define MATRIXPROPERTIES

#include <dumux/material/property_baseclasses.hh>

/**
 * \ingroup properties
 * \author Jochen Fritz
 */

namespace Dune
{

template<class G, class RT>
class Homogeneoussoil:Matrix2p<G,RT>
{
public:
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename G::ctype DT;
	enum {n=G::dimension, m=1};
	
	virtual const FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
	{
		return K_;
	}
	virtual double porosity(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
	{
		return 1.;
	}

	virtual double Sr_w(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T) const
	{
		return 0;
	}
	
	virtual double Sr_n(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T) const
	{
		return 0;
	}

	/* ATTENTION: define heat capacity per cubic meter! Be sure, that it corresponds to porosity!
			 * Best thing will be to define heatCap = (specific heatCapacity of material) * density * porosity*/
	virtual double heatCap(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
	{
		return 	790 /* spec. heat cap. of granite */
						* 2700 /* density of granite */
						* porosity(x, e, xi);
	}
	
	virtual double heatCond(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double sat) const
	{
		static const double lWater = 0.6;
		static const double lGranite = 2.8;
		double poro = porosity(x, e, xi);
		double lsat = pow(lGranite, (1-poro)) * pow(lWater, poro);
		double ldry = pow(lGranite, (1-poro));
		return ldry + sqrt(sat) * (ldry - lsat);
	}
		
	virtual std::vector<double> paramRelPerm(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T) const
	{
		// example for Brooks-Corey parameters
		std::vector<double> param(2);
		param[0] = 2.; // lambda
		param[1] = 0.; // entry-pressure
		return param;
	}
	
	virtual int relPermFlag(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
	{
		return 0;
	}
	
	Homogeneoussoil():Matrix2p<G,RT>()
	{
		for(int i = 0; i < n; i++)
			K_[i][i] = 1e-10;
	}
	
	~Homogeneoussoil()
	{}
	
private:
	FieldMatrix<DT,n,n> K_;
		
};

} // end namespace
#endif /*MATRIXPROPERTIES*/
