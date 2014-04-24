/*
  Copyright (C) 2014 by Andreas Lauser

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
*/
/*!
 * \file
 *
 * \copydoc Ewoms::EclipseWriter
 */
#ifndef EWOMS_ECLIPSE_WRITER_HH
#define EWOMS_ECLIPSE_WRITER_HH

#include "baseoutputwriter.hh"
#include "ertwrappers.hh"

#include <ewoms/disc/ecfv/ecfvdiscretization.hh>

#include <opm/material/Valgrind.hpp>

#if HAVE_DUNE_CORNERPOINT
#include <dune/grid/CpGrid.hpp>
#endif

#include <boost/algorithm/string.hpp>

#include <list>
#include <string>
#include <limits>
#include <sstream>
#include <fstream>
#include <type_traits>

// forward declaration of Dune::CpGrid
namespace Dune {
class CpGrid;
}

namespace Ewoms {
template <class TypeTag>
class EclipseWriter;

// we need to use specialization here because we need to call
// Dune::CpGrid specific methods (i.e. methods which are not part of
// the Dune Grid Interface and which never will be). This would cause
// the compiler to bail out if a grid different than Dune::CpGrid was
// used. (one could also take advantage of the HAVE_DUNE_CORNERPOINT
// macro, but even then the compiler would fail if dune-cornerpoint
// was available but a different grid type was chosen...)

// this is the stub class for non-cornerpoint grids
template <class TypeTag, class GridType>
class EclipseWriterHelper
{
    friend class EclipseWriter<TypeTag>;

    static void writeHeaders_(EclipseWriter<TypeTag> &writer)
    {
        OPM_THROW(std::logic_error,
                  "Eclipse binary output requires to use Dune::CpGrid");
    }
};

// this is the "real" code for cornerpoint grids
template <class TypeTag>
class EclipseWriterHelper<TypeTag, Dune::CpGrid>
{
    friend class EclipseWriter<TypeTag>;

    static void writeHeaders_(EclipseWriter<TypeTag> &writer)
    {
        typedef typename GET_PROP_TYPE(TypeTag, Discretization) Discretization;
        if (!std::is_same<Discretization, Ewoms::EcfvDiscretization<TypeTag> >::value)
            OPM_THROW(std::logic_error,
                      "Eclipse binary output only works for the element centered "
                      "finite volume discretization.");

#if ! HAVE_ERT
        OPM_THROW(std::logic_error,
                  "Eclipse binary output requires the ERT libraries");
#else

        // set the index of the first time step written to 0...
        writer.reportStepIdx_ = 0;

        char* egridRawFileName = ecl_util_alloc_filename(/*outputDir=*/"./",
                                                         writer.caseName().c_str(),
                                                         ECL_EGRID_FILE,
                                                         /*formatted=*/false, // -> write binary output
                                                         writer.reportStepIdx_);
        std::string egridFileName(egridRawFileName);
        std::free(egridRawFileName);

        ErtGrid ertGrid(writer.simulator_.gridManager().eclipseGrid());
        ertGrid.write(egridFileName, writer.reportStepIdx_);
#endif
    }
};

/*!
 * \brief Implements writing Eclipse binary output files.
 *
 * Caveats:
 * - For this class to do do anything meaningful, you need to have the
 *   ERT libraries with development headers installed and the ERT
 *   build system test must pass sucessfully.
 * - The only DUNE grid which is currently supported is Dune::CpGrid
 *   from the OPM module "dune-cornerpoint". Using another grid won't
 *   fail at compile time but you will provoke a fatal exception as
 *   soon as you try to write an Eclipse output file.
 * - This class requires to use the black oil model with the element
 *   centered finite volume discretization.
 * - MPI-parallel computations are not (yet?) supported.
 */
template <class TypeTag>
class EclipseWriter : public BaseOutputWriter
{
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, GridManager) GridManager;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum { dim = GridView::dimension };

    typedef BaseOutputWriter::Scalar Scalar;
    typedef BaseOutputWriter::Vector Vector;
    typedef BaseOutputWriter::ScalarBuffer ScalarBuffer;
    typedef BaseOutputWriter::VectorBuffer VectorBuffer;

    friend class EclipseWriterHelper<TypeTag, Grid>;

public:
    EclipseWriter(const Simulator &simulator)
        : simulator_(simulator)
        , gridView_(simulator_.gridView())
        , elementMapper_(gridView_)
        , vertexMapper_(gridView_)
    {
        reportStepIdx_ = 0;

        if (enableEclipseOutput_())
            EclipseWriterHelper<TypeTag, Grid>::writeHeaders_(*this);
    }

    ~EclipseWriter()
    { }

    /*!
     * \brief Returns the name of the simulation.
     *
     * This is the prefix of the files written to disk.
     */
    std::string caseName() const
    { return boost::to_upper_copy(simulator_.problem().name()); }

    /*!
     * \brief Updates the internal data structures after mesh
     *        refinement.
     *
     * If the grid changes between two calls of beginWrite(), this
     * method _must_ be called before the second beginWrite()!
     */
    void gridChanged()
    {
        elementMapper_.update();
        vertexMapper_.update();
    }

    /*!
     * \brief Called when ever a new time step or a new grid must be
     *        written.
     */
    void beginWrite(double t)
    {}

    /*
     * \brief Add a vertex-centered scalar field to the output.
     *
     * For the EclipseWriter, this method is a no-op which throws a
     * std::logic_error exception
     */
    void attachScalarVertexData(ScalarBuffer &buf, std::string name)
    {
        OPM_THROW(std::logic_error,
                  "The EclipseWriter can only write element based quantities!");
    }

    /*
     * \brief Add a vertex-centered vector field to the output.
     *
     * For the EclipseWriter, this method is a no-op which throws a
     * std::logic_error exception
     */
    void attachVectorVertexData(VectorBuffer &buf, std::string name)
    {
        OPM_THROW(std::logic_error,
                  "The EclipseWriter can only write element based quantities!");
    }

    /*!
     * \brief Add a scalar quantity to the output.
     *
     * The buffer must exist at least until the call to endWrite()
     * finishes. Modifying the buffer between the call to this method
     * and endWrite() results in _undefined behavior_.
     */
    void attachScalarElementData(ScalarBuffer &buf, std::string name)
    {
        attachedBuffers_.push_back(&buf);
    }

    /*
     * \brief Add a element-centered vector field to the output.
     *
     * For the EclipseWriter, this method is a no-op which throws a
     * std::logic_error exception
     */
    void attachVectorElementData(VectorBuffer &buf, std::string name)
    {
        OPM_THROW(std::logic_error,
                  "Currently, the EclipseWriter can only write scalar quantities!");
    }

    /*!
     * \brief Finalizes the current writer.
     *
     * This means that everything will be written to disk, except if
     * the onlyDiscard argument is true. In this case, no output is
     * written but the 'janitorial' jobs at the end of a time step are
     * still done.
     */
    void endWrite(bool onlyDiscard = false)
    {
        if (onlyDiscard) {
            // detach all buffers
            attachedBuffers_.resize(0);

            --reportStepIdx_;

            return;
        }

#warning "TODO: write the restart file!"

        // detach all buffers
        attachedBuffers_.resize(0);
    }

    /*!
     * \brief Write the multi-writer's state to a restart file.
     */
    template <class Restarter>
    void serialize(Restarter &res)
    {
        res.serializeSectionBegin("EclipseWriter");
        res.serializeStream() << reportStepIdx_ << "\n";
        res.serializeSectionEnd();
    }

    /*!
     * \brief Read the multi-writer's state from a restart file.
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        res.deserializeSectionBegin("EclipseWriter");
        res.deserializeStream() >> reportStepIdx_;
        res.deserializeSectionEnd();
    }

private:
    static bool enableEclipseOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableEclipseOutput); }

    // make sure the field is well defined if running under valgrind
    // and make sure that all values can be displayed by paraview
    void sanitizeBuffer_(std::vector<float> &b)
    {
        static bool warningPrinted = false;
        for (size_t i = 0; i < b.size(); ++i) {
            Valgrind::CheckDefined(b[i]);

            if (!warningPrinted && !std::isfinite(b[i])) {
                std::cerr << "WARNING: data field written to disk contains non-finite entries!\n";
                warningPrinted = true;
            }

            // set values which are too small to 0 to avoid possible
            // problems
            if (std::abs(b[i]) < std::numeric_limits<float>::min()) {
                b[i] = 0.0;
            }
        }
    }

    const Simulator& simulator_;
    const GridView gridView_;

    ElementMapper elementMapper_;
    VertexMapper vertexMapper_;

    double curTime_;
    int reportStepIdx_;

    std::list<ScalarBuffer*> attachedBuffers_;
};
} // namespace Ewoms

#endif
