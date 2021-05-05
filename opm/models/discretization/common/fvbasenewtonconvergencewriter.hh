// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Opm::FvBaseNewtonConvergenceWriter
 */
#ifndef EWOMS_FV_BASE_NEWTON_CONVERGENCE_WRITER_HH
#define EWOMS_FV_BASE_NEWTON_CONVERGENCE_WRITER_HH

#include <opm/models/io/vtkmultiwriter.hh>
#include <opm/models/utils/propertysystem.hh>

#include <iostream>

//! \cond SKIP_THIS
namespace Opm::Properties {

// forward declaration of the required property tags
template<class TypeTag, class MyTypeTag>
struct SolutionVector;
template<class TypeTag, class MyTypeTag>
struct GlobalEqVector;
template<class TypeTag, class MyTypeTag>
struct NewtonMethod;
template<class TypeTag, class MyTypeTag>
struct VtkOutputFormat;

} // namespace Opm::Properties
//! \endcond

namespace Opm {
/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief Writes the intermediate solutions during the Newton scheme
 *        for models using a finite volume discretization
 */
template <class TypeTag>
class FvBaseNewtonConvergenceWriter
{
    using GridView = GetPropType<TypeTag, Properties::GridView>;

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using NewtonMethod = GetPropType<TypeTag, Properties::NewtonMethod>;

    static const int vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;

public:
    FvBaseNewtonConvergenceWriter(NewtonMethod& nm)
        : newtonMethod_(nm)
    {
        timeStepIdx_ = 0;
        iteration_ = 0;
        vtkMultiWriter_ = 0;
    }

    ~FvBaseNewtonConvergenceWriter()
    { delete vtkMultiWriter_; }

    /*!
     * \brief Called by the Newton method before the actual algorithm
     *        is started for any given timestep.
     */
    void beginTimeStep()
    {
        ++timeStepIdx_;
        iteration_ = 0;
    }

    /*!
     * \brief Called by the Newton method before an iteration of the
     *        Newton algorithm is started.
     */
    void beginIteration()
    {
        ++ iteration_;
        if (!vtkMultiWriter_)
            vtkMultiWriter_ =
                new VtkMultiWriter(/*async=*/false,
                                   newtonMethod_.problem().gridView(),
                                   newtonMethod_.problem().outputDir(),
                                   "convergence");
        vtkMultiWriter_->beginWrite(timeStepIdx_ + iteration_ / 100.0);
    }

    /*!
     * \brief Write the Newton update to disk.
     *
     * Called after the linear solution is found for an iteration.
     *
     * \param uLastIter The solution vector of the previous iteration.
     * \param deltaU The negative difference between the solution
     *        vectors of the previous and the current iteration.
     */
    void writeFields(const SolutionVector& uLastIter,
                     const GlobalEqVector& deltaU)
    {
        try {
            newtonMethod_.problem().model().addConvergenceVtkFields(*vtkMultiWriter_,
                                                                    uLastIter,
                                                                    deltaU);
        }
        catch (...) {
            std::cout << "Oops: exception thrown on rank "
                      << newtonMethod_.problem().gridView().comm().rank()
                      << " while writing the convergence\n"  << std::flush;
        }
    }

    /*!
     * \brief Called by the Newton method after an iteration of the
     *        Newton algorithm has been completed.
     */
    void endIteration()
    { vtkMultiWriter_->endWrite(); }

    /*!
     * \brief Called by the Newton method after Newton algorithm
     *        has been completed for any given timestep.
     *
     * This method is called regardless of whether the Newton method
     * converged or not.
     */
    void endTimeStep()
    { iteration_ = 0; }

private:
    int timeStepIdx_;
    int iteration_;
    VtkMultiWriter *vtkMultiWriter_;
    NewtonMethod& newtonMethod_;
};

} // namespace Opm

#endif
