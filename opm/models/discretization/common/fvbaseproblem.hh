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
 * \copydoc Opm::FvBaseProblem
 */
#ifndef EWOMS_FV_BASE_PROBLEM_HH
#define EWOMS_FV_BASE_PROBLEM_HH

#include <dune/common/fvector.hh>

#include <opm/models/discretization/common/fvbaseparameters.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/discretization/common/restrictprolong.hh>

#include <opm/models/io/vtkmultiwriter.hh>
#include <opm/models/io/restart.hpp>

#include <opm/models/utils/simulatorutils.hpp>

#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>

namespace Opm::Properties {

template <class TypeTag, class MyTypeTag>
struct NewtonMethod;

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief Base class for all problems which use a finite volume spatial discretization.
 *
 * \note All quantities are specified assuming a threedimensional world. Problems
 *       discretized using 2D grids are assumed to be extruded by \f$1 m\f$ and 1D grids
 *       are assumed to have a cross section of \f$1m \times 1m\f$.
 */
template<class TypeTag>
class FvBaseProblem
{
private:
    using Implementation = GetPropType<TypeTag, Properties::Problem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

    static constexpr auto vtkOutputFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkOutputFormat>;

    using Model = GetPropType<TypeTag, Properties::Model>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using ThreadManager = GetPropType<TypeTag, Properties::ThreadManager>;
    using NewtonMethod = GetPropType<TypeTag, Properties::NewtonMethod>;

    using VertexMapper = GetPropType<TypeTag, Properties::VertexMapper>;
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;

    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Constraints = GetPropType<TypeTag, Properties::Constraints>;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;
    using VertexIterator = typename GridView::template Codim<dim>::Iterator;

    using CoordScalar = typename GridView::Grid::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

public:
    // the default restriction and prolongation for adaptation is simply an empty one
    using RestrictProlongOperator = EmptyRestrictProlong ;

private:
    // copying a problem is not a good idea
    FvBaseProblem(const FvBaseProblem&) = delete;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     *
     * \param simulator The time manager of the simulation
     * \param gridView The view on the DUNE grid which ought to be
     *                 used (normally the leaf grid view)
     */
    explicit FvBaseProblem(Simulator& simulator)
        : nextTimeStepSize_(0.0)
        , gridView_(simulator.gridView())
        , elementMapper_(gridView_, Dune::mcmgElementLayout())
        , vertexMapper_(gridView_, Dune::mcmgVertexLayout())
        , boundingBoxMin_(std::numeric_limits<double>::max())
        , boundingBoxMax_(-std::numeric_limits<double>::max())
        , simulator_(simulator)
    {
        // calculate the bounding box of the local partition of the grid view
        for (const auto& vertex : vertices(gridView_)) {
            for (unsigned i = 0; i < dim; ++i) {
                boundingBoxMin_[i] = std::min(boundingBoxMin_[i], vertex.geometry().corner(0)[i]);
                boundingBoxMax_[i] = std::max(boundingBoxMax_[i], vertex.geometry().corner(0)[i]);
            }
        }

        // communicate to get the bounding box of the whole domain
        for (unsigned i = 0; i < dim; ++i) {
            boundingBoxMin_[i] = gridView_.comm().min(boundingBoxMin_[i]);
            boundingBoxMax_[i] = gridView_.comm().max(boundingBoxMax_[i]);
        }

        if (enableVtkOutput_()) {
            const bool asyncVtkOutput =
                simulator_.gridView().comm().size() == 1 &&
                Parameters::Get<Parameters::EnableAsyncVtkOutput>();

            // asynchonous VTK output currently does not work in conjunction with grid
            // adaptivity because the async-IO code assumes that the grid stays
            // constant. complain about that case.
            const bool enableGridAdaptation = Parameters::Get<Parameters::EnableGridAdaptation>();
            if (asyncVtkOutput && enableGridAdaptation) {
                throw std::runtime_error("Asynchronous VTK output currently cannot be used "
                                         "at the same time as grid adaptivity");
            }

            defaultVtkWriter_ = std::make_unique<VtkMultiWriter>(asyncVtkOutput,
                                                                 gridView_,
                                                                 asImp_().outputDir(),
                                                                 asImp_().name());
        }
    }

    /*!
     * \brief Registers all available parameters for the problem and
     *        the model.
     */
    static void registerParameters()
    {
        Model::registerParameters();
        Parameters::Register<Parameters::MaxTimeStepSize<Scalar>>
            ("The maximum size to which all time steps are limited to [s]");
        Parameters::Register<Parameters::MinTimeStepSize<Scalar>>
            ("The minimum size to which all time steps are limited to [s]");
        Parameters::Register<Parameters::MaxTimeStepDivisions>
            ("The maximum number of divisions by two of the timestep size "
             "before the simulation bails out");
        Parameters::Register<Parameters::EnableAsyncVtkOutput>
            ("Dispatch a separate thread to write the VTK output");
        Parameters::Register<Parameters::ContinueOnConvergenceError>
            ("Continue with a non-converged solution instead of giving up "
             "if we encounter a time step size smaller than the minimum time "
             "step size.");
    }

    /*!
     * \brief Return if the storage term of the first iteration is identical to the storage
     *        term for the solution of the previous time step.
     *
     * This is only relevant if the storage cache is enabled and is usually the case,
     * i.e., this method only needs to be overwritten in rare corner cases.
     */
    bool recycleFirstIterationStorage() const
    { return true; }

    /*!
     * \brief Determine the directory for simulation output.
     *
     * The actual problem may chose to transform the value of the OutputDir parameter and
     * it can e.g. choose to create the directory on demand if it does not exist. The
     * default behaviour is to just return the OutputDir parameter and to throw an
     * exception if no directory with this name exists.
     */
    std::string outputDir() const
    {
        return simulatorOutputDir();
    }

    /*!
     * \brief Returns the string that is printed before the list of command line
     *        parameters in the help message.
     *
     * If the returned string is empty, no help message will be generated.
     */
    static std::string helpPreamble(int, const char **argv)
    {
        std::string desc = Implementation::briefDescription();
        if (!desc.empty()) {
            desc = desc + "\n";
        }

        return "Usage: " + std::string(argv[0]) + " [OPTIONS]\n" + desc;
    }

    /*!
     * \brief Returns a human readable description of the problem for the help message
     *
     * The problem description is printed as part of the --help message. It is optional
     * and should not exceed one or two lines of text.
     */
    static std::string briefDescription()
    { return ""; }

    // TODO (?): detailedDescription()

    /*!
     * \brief Handles positional command line parameters.
     *
     * Positional parameters are parameters that are not prefixed by any parameter name.
     *
     * \param seenParams The parameters which have already been seen in the current context
     * \param errorMsg If the positional argument cannot be handled, this is the reason why
     * \param argc The total number of command line parameters
     * \param argv The string value of the command line parameters
     * \param paramIdx The index of the positional parameter in the array of all parameters
     * \param posParamIdx The number of the positional parameter encountered so far
     *
     * \return The number of array entries which ought to be skipped before processing
     *         the next regular parameter. If this is less than 1, it indicated that the
     *         positional parameter was invalid.
     */
    static int handlePositionalParameter(std::function<void(const std::string&,
                                                            const std::string&)>,
                                         std::set<std::string>&,
                                         std::string& errorMsg,
                                         int,
                                         const char** argv,
                                         int paramIdx,
                                         int)
    {
        errorMsg = std::string("Illegal parameter \"") + argv[paramIdx] + "\".";
        return 0;
    }

    /*!
     * \brief Called by the Opm::Simulator in order to initialize the problem.
     *
     * If you overload this method don't forget to call ParentType::finishInit()
     */
    void finishInit()
    {}

    /*!
     * \brief Allows to improve the performance by prefetching all data which is
     *        associated with a given element.
     */
    void prefetch(const Element&) const
    {
        // do nothing by default
    }

    /*!
     * \brief Handle changes of the grid
     */
    void gridChanged()
    {
        elementMapper_.update(gridView_);
        vertexMapper_.update(gridView_);

        if (enableVtkOutput_()) {
            defaultVtkWriter_->gridChanged();
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a boundary segment.
     *
     * \param values Stores the fluxes over the boundary segment.
     * \param context The object representing the execution context from
     *                which this method is called.
     * \param spaceIdx The local index of the spatial entity which represents the boundary segment.
     * \param timeIdx The index used for the time discretization
     */
    template <class Context>
    void boundary(BoundaryRateVector&,
                  const Context&,
                  unsigned,
                  unsigned) const
    { throw std::logic_error("Problem does not provide a boundary() method"); }

    /*!
     * \brief Evaluate the constraints for a control volume.
     *
     * \param constraints Stores the values of the primary variables at a
     *                    given spatial and temporal location.
     * \param context The object representing the execution context from
     *                which this method is called.
     * \param spaceIdx The local index of the spatial entity which represents the boundary segment.
     * \param timeIdx The index used for the time discretization
     */
    template <class Context>
    void constraints(Constraints&,
                     const Context&,
                     unsigned,
                     unsigned) const
    { throw std::logic_error("Problem does not provide a constraints() method"); }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param rate Stores the values of the volumetric creation/anihilition
     *             rates of the conserved quantities.
     * \param context The object representing the execution context from which
     *                this method is called.
     * \param spaceIdx The local index of the spatial entity which represents
     *                 the boundary segment.
     * \param timeIdx The index used for the time discretization
     */
    template <class Context>
    void source(RateVector&,
                const Context&,
                unsigned,
                unsigned) const
    { throw std::logic_error("Problem does not provide a source() method"); }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values Stores the primary variables.
     * \param context The object representing the execution context from which
     *                this method is called.
     * \param spaceIdx The local index of the spatial entity which represents
     *                 the boundary segment.
     * \param timeIdx The index used for the time discretization
     */
    template <class Context>
    void initial(PrimaryVariables&,
                 const Context&,
                 unsigned,
                 unsigned) const
    { throw std::logic_error("Problem does not provide a initial() method"); }

    /*!
     * \brief Return how much the domain is extruded at a given sub-control volume.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     *
     * \param context The object representing the execution context from which
     *                this method is called.
     * \param spaceIdx The local index of the spatial entity which represents
     *                 the boundary segment.
     * \param timeIdx The index used for the time discretization
     */
    template <class Context>
    Scalar extrusionFactor(const Context&,
                           unsigned,
                           unsigned) const
    { return asImp_().extrusionFactor(); }

    Scalar extrusionFactor() const
    { return 1.0; }

    /*!
     * \brief Callback used by the model to indicate that the initial solution has been
     *        determined for all degrees of freedom.
     */
    void initialSolutionApplied()
    {}

    /*!
     * \brief Called at the beginning of an simulation episode.
     */
    void beginEpisode()
    {}

    /*!
     * \brief Called by the simulator before each time integration.
     */
    void beginTimeStep()
    {}

    /*!
     * \brief Called by the simulator before each Newton-Raphson iteration.
     */
    void beginIteration()
    {}

    /*!
     * \brief Called by the simulator after each Newton-Raphson update.
     */
    void endIteration()
    {}

    /*!
     * \brief Called by the simulator after each time integration.
     *
     * This method is intended to do some post processing of the
     * solution. (e.g., some additional output)
     */
    void endTimeStep()
    {}

    /*!
     * \brief Called when the end of an simulation episode is reached.
     *
     * Typically, a new episode is started in this method.
     */
    void endEpisode()
    {
        std::cerr << "The end of episode " << simulator().episodeIndex() + 1 << " has been "
                  << "reached, but the problem does not override the endEpisode() method. "
                  << "Doing nothing!\n";
    }

    /*!
     * \brief Called after the simulation has been run sucessfully.
     */
    void finalize()
    {
        const auto& executionTimer = simulator().executionTimer();

        const Scalar executionTime = executionTimer.realTimeElapsed();
        const Scalar setupTime = simulator().setupTimer().realTimeElapsed();
        const Scalar prePostProcessTime = simulator().prePostProcessTimer().realTimeElapsed();
        const Scalar localCpuTime = executionTimer.cpuTimeElapsed();
        const Scalar globalCpuTime = executionTimer.globalCpuTimeElapsed();
        const Scalar writeTime = simulator().writeTimer().realTimeElapsed();
        const Scalar linearizeTime = simulator().linearizeTimer().realTimeElapsed();
        const Scalar solveTime = simulator().solveTimer().realTimeElapsed();
        const Scalar updateTime = simulator().updateTimer().realTimeElapsed();
        const unsigned numProcesses = static_cast<unsigned>(this->gridView().comm().size());
        const unsigned threadsPerProcess = ThreadManager::maxThreads();
        if (gridView().comm().rank() == 0) {
            std::cout << std::setprecision(3)
                      << "Simulation of problem '" << asImp_().name() << "' finished.\n"
                      << "\n"
                      << "------------------------ Timing ------------------------\n"
                      << "Setup time: " << setupTime << " seconds"
                      << humanReadableTime(setupTime)
                      << ", " << setupTime / (executionTime + setupTime) * 100 << "%\n"
                      << "Simulation time: " << executionTime << " seconds"
                      << humanReadableTime(executionTime)
                      << ", " << executionTime / (executionTime + setupTime) * 100 << "%\n"
                      << "    Linearization time: " << linearizeTime << " seconds"
                      << humanReadableTime(linearizeTime)
                      << ", " << linearizeTime / executionTime * 100 << "%\n"
                      << "    Linear solve time: "  << solveTime << " seconds"
                      << humanReadableTime(solveTime)
                      << ", " << solveTime / executionTime * 100 << "%\n"
                      << "    Newton update time: "  << updateTime << " seconds"
                      << humanReadableTime(updateTime)
                      << ", " << updateTime / executionTime * 100 << "%\n"
                      << "    Pre/postprocess time: "  << prePostProcessTime << " seconds"
                      << humanReadableTime(prePostProcessTime)
                      << ", " << prePostProcessTime / executionTime * 100 << "%\n"
                      << "    Output write time: "  << writeTime << " seconds"
                      << humanReadableTime(writeTime)
                      << ", " << writeTime / executionTime * 100 << "%\n"
                      << "First process' simulation CPU time: "  << localCpuTime << " seconds"
                      <<  humanReadableTime(localCpuTime) << "\n"
                      << "Number of processes: " << numProcesses << "\n"
                      << "Threads per processes: " << threadsPerProcess << "\n"
                      << "Total CPU time: " << globalCpuTime << " seconds"
                      << humanReadableTime(globalCpuTime) << "\n"
                      << "\n"
                      << "----------------------------------------------------------------\n"
                      << std::endl;
        }
    }

    /*!
     * \brief Called by Opm::Simulator in order to do a time
     *        integration on the model.
     */
    void timeIntegration()
    {
        const unsigned maxFails = asImp_().maxTimeIntegrationFailures();
        Scalar minTimeStep = asImp_().minTimeStepSize();

        std::string errorMessage;
        for (unsigned i = 0; i < maxFails; ++i) {
            bool converged = model().update();
            if (converged) {
                return;
            }

            const Scalar dt = simulator().timeStepSize();
            Scalar nextDt = dt / 2.0;
            if (dt < minTimeStep * (1.0 + 1e-9)) {
                if (asImp_().continueOnConvergenceError()) {
                    if (gridView().comm().rank() == 0) {
                        std::cout << "Newton solver did not converge with minimum time step of "
                                  << dt << " seconds. Continuing with unconverged solution!\n"
                                  << std::flush;
                    }
                    return;
                }
                else {
                    errorMessage = "Time integration did not succeed with the minumum time step size of " +
                                   std::to_string(double(minTimeStep)) + " seconds. Giving up!";
                    break; // give up: we can't make the time step smaller anymore!
                }
            }
            else if (nextDt < minTimeStep) {
                nextDt = minTimeStep;
            }
            simulator().setTimeStepSize(nextDt);

            // update failed
            if (gridView().comm().rank() == 0) {
                std::cout << "Newton solver did not converge with "
                          << "dt=" << dt << " seconds. Retrying with time step of "
                          << nextDt << " seconds\n" << std::flush;
            }
        }

        if (errorMessage.empty()) {
            errorMessage = "Newton solver didn't converge after " +
                           std::to_string(maxFails) + " time-step divisions. dt=" +
                           std::to_string(double(simulator().timeStepSize()));
        }
        throw std::runtime_error(errorMessage);
    }

    /*!
     * \brief Returns the minimum allowable size of a time step.
     */
    Scalar minTimeStepSize() const
    { return Parameters::Get<Parameters::MinTimeStepSize<Scalar>>(); }

    /*!
     * \brief Returns the maximum number of subsequent failures for the time integration
     *        before giving up.
     */
    unsigned maxTimeIntegrationFailures() const
    { return Parameters::Get<Parameters::MaxTimeStepDivisions>(); }

    /*!
     * \brief Returns if we should continue with a non-converged solution instead of
     *        giving up if we encounter a time step size smaller than the minimum time
     *        step size.
     */
    bool continueOnConvergenceError() const
    { return Parameters::Get<Parameters::ContinueOnConvergenceError>(); }

    /*!
     * \brief Impose the next time step size to be used externally.
     */
    void setNextTimeStepSize(Scalar dt)
    { nextTimeStepSize_ = dt; }

    /*!
     * \brief Called by Opm::Simulator whenever a solution for a
     *        time step has been computed and the simulation time has
     *        been updated.
     */
    Scalar nextTimeStepSize() const
    {
        if (nextTimeStepSize_ > 0.0) {
            return nextTimeStepSize_;
        }

        Scalar dtNext = std::min(Parameters::Get<Parameters::MaxTimeStepSize<Scalar>>(),
                                 newtonMethod().suggestTimeStepSize(simulator().timeStepSize()));

        if (dtNext < simulator().maxTimeStepSize() &&
            simulator().maxTimeStepSize() < dtNext * 2)
        {
            dtNext = simulator().maxTimeStepSize() / 2 * 1.01;
        }

        return dtNext;
    }

    /*!
     * \brief Returns true if a restart file should be written to
     *        disk.
     *
     * The default behavior is to write one restart file every 10 time
     * steps. This method should be overwritten by the
     * implementation if the default behavior is deemed insufficient.
     */
    bool shouldWriteRestartFile() const
    {
        return simulator().timeStepIndex() > 0 &&
               (simulator().timeStepIndex() % 10 == 0);
    }

    /*!
     * \brief Returns true if the current solution should be written to
     *        disk (i.e. as a VTK file)
     *
     * The default behavior is to write out the solution for every
     * time step. This method is should be overwritten by the
     * implementation if the default behavior is deemed insufficient.
     */
    bool shouldWriteOutput() const
    { return true; }

    /*!
     * \brief Called by the simulator after everything which can be
     *        done about the current time step is finished and the
     *        model should be prepared to do the next time integration.
     */
    void advanceTimeLevel()
    { model().advanceTimeLevel(); }

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     * It is highly recommend to overwrite this method in the concrete
     * problem which is simulated.
     */
    std::string name() const
    { return "sim"; }

    /*!
     * \brief The GridView which used by the problem.
     */
    const GridView& gridView() const
    { return gridView_; }

    /*!
     * \brief The coordinate of the corner of the GridView's bounding
     *        box with the smallest values.
     */
    const GlobalPosition& boundingBoxMin() const
    { return boundingBoxMin_; }

    /*!
     * \brief The coordinate of the corner of the GridView's bounding
     *        box with the largest values.
     */
    const GlobalPosition& boundingBoxMax() const
    { return boundingBoxMax_; }

    /*!
     * \brief Returns the mapper for vertices to indices.
     */
    const VertexMapper& vertexMapper() const
    { return vertexMapper_; }

    /*!
     * \brief Returns the mapper for elements to indices.
     */
    const ElementMapper& elementMapper() const
    { return elementMapper_; }

    /*!
     * \brief Returns Simulator object used by the simulation
     */
    Simulator& simulator()
    { return simulator_; }

    /*!
     * \copydoc simulator()
     */
    const Simulator& simulator() const
    { return simulator_; }

    /*!
     * \brief Returns numerical model used for the problem.
     */
    Model& model()
    { return simulator_.model(); }

    /*!
     * \copydoc model()
     */
    const Model& model() const
    { return simulator_.model(); }

    /*!
     * \brief Returns object which implements the Newton method.
     */
    NewtonMethod& newtonMethod()
    { return model().newtonMethod(); }

    /*!
     * \brief Returns object which implements the Newton method.
     */
    const NewtonMethod& newtonMethod() const
    { return model().newtonMethod(); }
    // \}

    /*!
     * \brief return restriction and prolongation operator
     * \note This method has to be overloaded by the implementation.
     */
    RestrictProlongOperator restrictProlongOperator()
    {
        return RestrictProlongOperator();
    }

    /*!
     * \brief Mark grid cells for refinement or coarsening
     * \note This method has to be overloaded in derived classes to proper implement
     *       marking of grid entities.
     *
     * \return number of marked cells (default is 0)
     */
    unsigned markForGridAdaptation()
    {
        return 0;
    }

    /*!
     * \brief This method writes the complete state of the problem
     *        to the harddisk.
     *
     * The file will start with the prefix returned by the name()
     * method, has the current time of the simulation clock in it's
     * name and uses the extension <tt>.ers</tt>. (Ewoms ReStart
     * file.)  See Opm::Restart for details.
     *
     * \tparam Restarter The serializer type
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void serialize(Restarter& res)
    {
        if (enableVtkOutput_()) {
            defaultVtkWriter_->serialize(res);
        }
    }

    /*!
     * \brief This method restores the complete state of the problem
     *        from disk.
     *
     * It is the inverse of the serialize() method.
     *
     * \tparam Restarter The deserializer type
     *
     * \param res The deserializer object
     */
    template <class Restarter>
    void deserialize(Restarter& res)
    {
        if (enableVtkOutput_()) {
            defaultVtkWriter_->deserialize(res);
        }
    }

    /*!
     * \brief Write the relevant secondary variables of the current
     *        solution into an VTK output file.
     *
     * \param verbose Specify if a message should be printed whenever a file is written
     */
    void writeOutput(bool verbose = true)
    {
        if (!enableVtkOutput_()) {
            return;
        }

        if (verbose && gridView().comm().rank() == 0) {
            std::cout << "Writing visualization results for the current time step.\n"
                      << std::flush;
        }

        // calculate the time _after_ the time was updated
        const Scalar t = simulator().time() + simulator().timeStepSize();

        defaultVtkWriter_->beginWrite(t);
        model().prepareOutputFields();
        model().appendOutputFields(*defaultVtkWriter_);
        defaultVtkWriter_->endWrite(false);
    }

    /*!
     * \brief Method to retrieve the VTK writer which should be used
     *        to write the default ouput after each time step to disk.
     */
    VtkMultiWriter& defaultVtkWriter() const
    { return *defaultVtkWriter_; }

protected:
    Scalar nextTimeStepSize_;

    bool enableVtkOutput_() const
    { return Parameters::Get<Parameters::EnableVtkOutput>(); }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    //! \copydoc asImp_()
    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

    // Grid management stuff
    const GridView gridView_;
    ElementMapper elementMapper_;
    VertexMapper vertexMapper_;
    GlobalPosition boundingBoxMin_;
    GlobalPosition boundingBoxMax_;

    // Attributes required for the actual simulation
    Simulator& simulator_;
    std::unique_ptr<VtkMultiWriter> defaultVtkWriter_{};
};

} // namespace Opm

#endif
