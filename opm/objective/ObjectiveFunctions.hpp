/*!
 * \file
 * \copydoc Opm::EclStone1Material
 */
#ifndef OPM_OBJECTIVEFUNCTIONS_HPP
#define  OPM_OBJECTIVEFUNCTIONS_HPP
#include <vector>
#include <cmath>
namespace Opm {
namespace Objectives {
class TotalWaterRate
{
    public:
    template<class RateType>
    static RateType objective(const std::vector<RateType>& rate,double t, double dt){
        assert(rate.size()>1);
        return rate[1]*dt;
    }
};
class OilProduction
{
    public:

    template<class RateType>
    RateType objective(const std::vector<RateType>& rate,double t, double dt){
        assert(rate.size()>1);
        assert(rate[0]>=0);
        auto prod = rate[0]*dt;
        double alph = t/tau_;
        double fac = std::pow(1+d_,-alph);
        return prod*fac;
    }
    private:
        static constexpr double d_ = 0.08;
        static constexpr double tau_ = 365*24*60*60;
};

class Olympus
{
    //http://www.isapp2.com/downloads/problem-statement.pdf
    public:

    template<class RateType>
    RateType objective(const std::vector<RateType>& rate,double t, double dt){
        assert(rate.size()==3);
        assert(rate[0]>=0);
        // oil production
        auto prod = rate[0]*dt;
        // water production
        if(rate[1]>0){
            prod -= rate[1]*rwp_;
        }else{
            assert(rate [1]<=0);
            prod += rate[1]*rwi_;
        }
        double alph = t/tau_;
        double fac = std::pow(1+d_,-alph);
        return prod*fac;
    }
    private:
        static constexpr double bbl_ = 0.1590;
        static constexpr double d_ = 0.08;
        static constexpr double tau_ = 365*24*60*60;
        static constexpr double rop_ = 45*bbl_;
        static constexpr double rwp_ = 6*bbl_;
        static constexpr double rwi_ = 2*bbl_;
};

}
}
#endif // OPM_OBJECTIVEFUNCTIONS_HPP
