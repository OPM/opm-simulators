/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <opm/core/fluid/SaturationPropsBasic.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <iostream>

namespace Opm
{


    // ---------- Helper functions ----------

    namespace {

	struct KrFunConstant
	{
	    double kr(double)
	    {
		return 1.0;
	    }
	    double dkrds(double)
	    {
		return 0.0;
	    }
	};

	struct KrFunLinear
	{
	    double kr(double s)
	    {
		return s;
	    }
	    double dkrds(double)
	    {
		return 1.0;
	    }
	};

	struct KrFunQuadratic
	{
	    double kr(double s)
	    {
		return s*s;
	    }
	    double dkrds(double s)
	    {
		return 2.0*s;
	    }
	};


	template <class Fun>
	static inline void evalAllKrDeriv(const int n, const int np,
					  const double* s, double* kr, double* dkrds, Fun fun)
	{
	    if (dkrds == 0) {
// #pragma omp parallel for
		for (int i = 0; i < n*np; ++i) {
		    kr[i] = fun.kr(s[i]);
		}
		return;
	    }
// #pragma omp parallel for
	    for (int i = 0; i < n; ++i) {
		std::fill(dkrds + i*np*np, dkrds + (i+1)*np*np, 0.0);
		for (int phase = 0; phase < np; ++phase) {
		    kr[i*np + phase] = fun.kr(s[i*np + phase]);
		    // Only diagonal elements in derivative.
		    dkrds[i*np*np + phase*np + phase] = fun.dkrds(s[i*np + phase]);
		}
	    }
	}


    } // anon namespace



    // ---------- Class methods ----------



    /// Default constructor.
    SaturationPropsBasic::SaturationPropsBasic()
    {
    }




    /// Initialize from parameters.
    void SaturationPropsBasic::init(const parameter::ParameterGroup& param)
    {
	int num_phases = param.getDefault("num_phases", 2);
	if (num_phases > 2 || num_phases < 1) {
	    THROW("SaturationPropsBasic::init() illegal num_phases: " << num_phases);
	}
        num_phases_ = num_phases;
	std::string rpf = param.getDefault("relperm_func", std::string("Linear"));
	if (rpf == "Constant") {
	    relperm_func_ = Constant;
	} else if (rpf == "Linear") {
	    relperm_func_ = Linear;
	} else if (rpf == "Quadratic") {
	    relperm_func_ = Quadratic;
	} else {
	    THROW("SaturationPropsBasic::init() illegal relperm_func: " << rpf);
	}
    }




    /// \return   P, the number of phases.
    int SaturationPropsBasic::numPhases() const
    {
	return num_phases_;
    }




    /// Relative permeability.
    /// \param[in]  n      Number of data points.
    /// \param[in]  s      Array of nP saturation values.
    /// \param[out] kr     Array of nP relperm values, array must be valid before calling.
    /// \param[out] dkrds  If non-null: array of nP^2 relperm derivative values,
    ///                    array must be valid before calling.
    ///                    The P^2 derivative matrix is
    ///                           m_{ij} = \frac{dkr_i}{ds^j},
    ///                    and is output in Fortran order (m_00 m_10 m_20 m01 ...)
    void SaturationPropsBasic::relperm(const int n,
				       const double* s,
				       double* kr,
				       double* dkrds) const
    {
	switch (relperm_func_) {
	case Constant:
	    {
		evalAllKrDeriv(n, num_phases_, s, kr, dkrds, KrFunConstant());
		break;
	    }
	case Linear:
	    {
		evalAllKrDeriv(n, num_phases_, s, kr, dkrds, KrFunLinear());
		break;
	    }
	case Quadratic:
	    {
		evalAllKrDeriv(n, num_phases_, s, kr, dkrds, KrFunQuadratic());
		break;
	    }
	default:
	    THROW("SaturationPropsBasic::relperm() unhandled relperm func type: " << relperm_func_);
	}
    }




    /// Capillary pressure.
    /// \param[in]  n      Number of data points.
    /// \param[in]  s      Array of nP saturation values.
    /// \param[out] pc     Array of nP capillary pressure values, array must be valid before calling.
    /// \param[out] dpcds  If non-null: array of nP^2 derivative values,
    ///                    array must be valid before calling.
    ///                    The P^2 derivative matrix is
    ///                           m_{ij} = \frac{dpc_i}{ds^j},
    ///                    and is output in Fortran order (m_00 m_10 m_20 m01 ...)
    void SaturationPropsBasic::capPress(const int n,
					const double* /*s*/,
					double* pc,
					double* dpcds) const
    {
	std::fill(pc, pc + num_phases_*n, 0.0);
        if (dpcds) {
	    std::fill(dpcds, dpcds + num_phases_*num_phases_*n, 0.0);
        }
    }



    /// Obtain the range of allowable saturation values.
    /// \param[in]  n      Number of data points.
    /// \param[out] smin   Array of nP minimum s values, array must be valid before calling.
    /// \param[out] smax   Array of nP maximum s values, array must be valid before calling.
    void SaturationPropsBasic::satRange(const int n,
					double* smin,
					double* smax) const
    {
	std::fill(smin, smin + num_phases_*n, 0.0);
	std::fill(smax, smax + num_phases_*n, 1.0);
    }



} // namespace Opm


