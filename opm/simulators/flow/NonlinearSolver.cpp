/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2015 Statoil ASA.

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
#include <config.h>
#include <opm/simulators/flow/NonlinearSolver.hpp>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <opm/common/ErrorMacros.hpp>

#include <cmath>
#include <stdexcept>

namespace Opm::detail {

template<class Scalar>
void detectOscillations(const std::vector<std::vector<Scalar>>& residualHistory,
                        const int it, const int numPhases, const Scalar relaxRelTol,
                        const int minimumOscillatingPhases,
                        bool& oscillate, bool& stagnate)
{
    // The detection of oscillation in two primary variable results in the report of the detection
    // of oscillation for the solver.
    // Only the saturations are used for oscillation detection for the black oil model.
    // Stagnate is not used for any treatment here.

    if (it < 2) {
        oscillate = false;
        stagnate = false;
        return;
    }

    stagnate = true;
    int oscillatePhase = 0;
    const auto& F0 = residualHistory[it];
    const auto& F1 = residualHistory[it - 1];
    const auto& F2 = residualHistory[it - 2];
    for (int p = 0; p < numPhases; ++p) {
        const Scalar d1 = std::abs((F0[p] - F2[p]) / F0[p]);
        const Scalar d2 = std::abs((F0[p] - F1[p]) / F0[p]);

        oscillatePhase += (d1 < relaxRelTol) && (relaxRelTol < d2);

        // Process is 'stagnate' unless at least one phase
        // exhibits significant residual change.
        stagnate = (stagnate && !(std::abs((F1[p] - F2[p]) / F2[p]) > 1.0e-3));
    }

    oscillate = (oscillatePhase >= minimumOscillatingPhases);
}

template <class BVector, class Scalar>
void stabilizeNonlinearUpdate(BVector& dx, BVector& dxOld,
                              const Scalar omega,
                              NonlinearRelaxType relaxType)
{
    // The dxOld is updated with dx.
    // If omega is equal to 1., no relaxtion will be appiled.

    BVector tempDxOld = dxOld;
    dxOld = dx;

    switch (relaxType) {
    case NonlinearRelaxType::Dampen: {
        if (omega == 1.) {
            return;
        }
        std::transform(dx.begin(), dx.end(), dx.begin(),
                       [omega](const auto d) { return d*omega; });
        return;
    }
    case NonlinearRelaxType::SOR: {
        if (omega == 1.) {
            return;
        }
        for (auto i = 0*dx.size(); i < dx.size(); ++i) {
            dx[i] *= omega;
            tempDxOld[i] *= (1.-omega);
            dx[i] += tempDxOld[i];
        }
        return;
    }
    default:
        OPM_THROW(std::runtime_error, "Can only handle Dampen and SOR relaxation type.");
    }

    return;
}

template<class Scalar, int Size>
using BV = Dune::BlockVector<Dune::FieldVector<Scalar,Size>>;

#define INSTANTIATE(T,Size)                                                 \
    template void stabilizeNonlinearUpdate(BV<T,Size>&, BV<T,Size>&,        \
                                           const T, NonlinearRelaxType);

#define INSTANTIATE_TYPE(T)                                                 \
    template void detectOscillations(const std::vector<std::vector<T>>&,    \
                                     const int, const int, const T,         \
                                     const int, bool&, bool&);              \
    INSTANTIATE(T,1)                                                        \
    INSTANTIATE(T,2)                                                        \
    INSTANTIATE(T,3)                                                        \
    INSTANTIATE(T,4)                                                        \
    INSTANTIATE(T,5)                                                        \
    INSTANTIATE(T,6)

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm::detail

namespace Opm {

template<class Scalar>
NonlinearSolverParameters<Scalar>::
NonlinearSolverParameters()
{
    // set default values
    reset();

    // overload with given parameters
    relaxMax_ = Parameters::Get<Parameters::NewtonMaxRelax<Scalar>>();

    const auto& relaxationTypeString = Parameters::Get<Parameters::NewtonRelaxationType>();
    if (relaxationTypeString == "dampen") {
        relaxType_ = NonlinearRelaxType::Dampen;
    } else if (relaxationTypeString == "sor") {
        relaxType_ = NonlinearRelaxType::SOR;
    } else {
        OPM_THROW(std::runtime_error,
                  "Unknown Relaxtion Type " + relaxationTypeString);
    }
}

template<class Scalar>
void NonlinearSolverParameters<Scalar>::
reset()
{
    // default values for the solver parameters
    relaxType_ = NonlinearRelaxType::Dampen;
    relaxMax_ = 0.5;
    relaxIncrement_ = 0.1;
    relaxRelTol_ = 0.2;
}

template<class Scalar>
void NonlinearSolverParameters<Scalar>::
registerParameters()
{
    Parameters::Register<Parameters::NewtonMaxRelax<Scalar>>
        ("The maximum relaxation factor of a Newton iteration");
    Parameters::Register<Parameters::NewtonRelaxationType>
        ("The type of relaxation used by Newton method");
}

template struct NonlinearSolverParameters<double>;

#if FLOW_INSTANTIATE_FLOAT
template struct NonlinearSolverParameters<float>;
#endif

}
