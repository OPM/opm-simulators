/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 IRIS AS

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

#ifndef OPM_NEWTONITERATIONBLACKOILSIMPLE_HEADER_INCLUDED
#define OPM_NEWTONITERATIONBLACKOILSIMPLE_HEADER_INCLUDED


#include <opm/autodiff/NewtonIterationBlackoilInterface.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <memory>

namespace Opm
{

    /// This class solves the fully implicit black-oil system by
    /// simply concatenating the Jacobian matrices and passing the
    /// resulting system to a linear solver. The linear solver used
    /// can be passed in as a constructor argument.
    class NewtonIterationBlackoilSimple : public NewtonIterationBlackoilInterface
    {
    public:
        /// Construct a system solver.
        /// \param[in] param   parameters controlling the behaviour and
        ///                    choice of linear solver.
        /// \param[in] parallelInformation In the case of a parallel run
        ///                    with dune-istl the information about the parallelization.
        NewtonIterationBlackoilSimple(const ParameterGroup& param,
                                      const boost::any& parallelInformation=boost::any());

        /// Solve the system of linear equations Ax = b, with A being the
        /// combined derivative matrix of the residual and b
        /// being the residual itself.
        /// \param[in] residual   residual object containing A and b.
        /// \return               the solution x
        virtual SolutionVector computeNewtonIncrement(const LinearisedBlackoilResidual& residual) const;

        /// \copydoc NewtonIterationBlackoilInterface::iterations
        virtual int iterations () const { return iterations_; }

        /// \copydoc NewtonIterationBlackoilInterface::parallelInformation
        virtual const boost::any& parallelInformation() const;

    private:
        std::unique_ptr<LinearSolverInterface> linsolver_;
        mutable int iterations_;
        boost::any parallelInformation_;
    };

} // namespace Opm


#endif // OPM_NEWTONITERATIONBLACKOILSIMPLE_HEADER_INCLUDED
