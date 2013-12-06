/*
Author : Liu Ming
Data   : 2013-11-28  in Oslo
Email  : miliu@statoil.com

Properties for incompressible immiscible two-phase flow
*/

#ifndef OPM_INCOMPROPSAD_HEADER_INCLUDED
#define  OPM_INCOMPROPSAD_HEADER_INCLUDED
#include <opm/core/props/IncompPropertiesBasic.hpp>
#include <opm/core/props/rock/RockBasic.hpp>
#include <opm/polymer/fullyimplicit/AutoDiffBlock.hpp>
#include <opm/core/props/pvt/PvtPropertiesBasic.hpp>
#include <opm/core/props/satfunc/SaturationPropsBasic.hpp>


namespace Opm
{

//    class BlackoilPhases;

    class IncompropsAd : public IncompPropertiesBasic
    {
    public:
        IncomPropsAd(const parameter::ParameterGroup& param,
                     const int dim,
                     const int num_cells);
        IncompPropsAd(const int num_phases,
                      const SaturationPropsBasic::RelPermFunc& relpermfunc,
                      const std::vector<double>& rho,
                      const std::vector<double>& mu,
                      const double porosity,
                      const double permeability,
                      const int dim,
                      const int num_cells);

        ~IncompPropsAd();

        typedef AutoDiffBlock<double> ADB;
        typedef ADB::V V;
        V relperm(const V& sw,
                  const V& so,
                  const std::vector<int>& cells);
        ADB relperm(const ADB& sw,
                    const ADB& so,
                    const std::vector<int>& cells);

        V capPress(const V& sw,
                   const V& so,
                   const std::vector<int>& cells);
        ADB capPress(const ADB& sw,
                     const ADB& so,
                     const std::vector<int>& cells);

    }
}

#endif//  OPM_INCOMPROPSAD_HEADER_INCLUDED
